"""
Diagnose slow or stalled ``import f4enix``.

The script never imports f4enix itself. Every measurement runs in a fresh
subprocess with a timeout, so a hang is caught and its stack is printed
(via ``faulthandler``) instead of freezing this script.

Usage (Linux or Windows, stdlib only)::

    python scripts/diagnose_import.py
    python scripts/diagnose_import.py --cwd /path/where/notebook/lives
    python scripts/diagnose_import.py --python C:/path/to/kernel/python.exe
    python scripts/diagnose_import.py --with-dead-display      # Linux only
    python scripts/diagnose_import.py --output report.txt

Run it from (or with ``--cwd`` pointing to) the folder where the import is
slow, e.g. the notebook folder on a network share or OneDrive, and with the
same interpreter the Jupyter kernel uses (``--python``).

What it reports
---------------
1. Environment: OS, interpreter, relevant env variables, filesystem type of
   the working directory and of every ``sys.path`` entry (network shares,
   OneDrive), matplotlib cache status, Windows Defender status (best effort).
2. A/B scenarios of ``import f4enix`` (timed, with timeout):
   baseline, repeated baseline (cold vs warm caches / antivirus), cwd removed
   from ``sys.path``, fresh matplotlib font cache, Agg backend, DISPLAY unset.
3. Per-dependency import times, each in its own subprocess.
4. Loaded OpenMP / BLAS runtimes (mixing libiomp and libomp can deadlock).
5. The 15 slowest modules (``-X importtime``) of the baseline run and the
   faulthandler stack of any run that timed out.

Diagnosing inside Jupyter
-------------------------
Put this in the first cell, before ``import f4enix``::

    import faulthandler, sys
    faulthandler.dump_traceback_later(15, repeat=True, file=sys.__stderr__)

While the cell hangs, the stack of the kernel is printed every 15 s to the
Jupyter server log (VS Code: Output panel -> "Jupyter"). Cancel with
``faulthandler.cancel_dump_traceback_later()``. Alternatively
``pip install py-spy`` and run ``py-spy dump --pid <kernel pid>`` while it
hangs. Compare ``%time import f4enix`` in the kernel with this script.
"""

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or -
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at:
    https://eupl.eu/1.2/en/
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""


import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field

IS_WIN = sys.platform.startswith("win")
IS_LINUX = sys.platform.startswith("linux")

ENV_VARS = [
    "DISPLAY",
    "WAYLAND_DISPLAY",
    "MPLBACKEND",
    "MPLCONFIGDIR",
    "XDG_CACHE_HOME",
    "PYTHONPATH",
    "PYTHONSAFEPATH",
    "PYTHONHOME",
    "CONDA_PREFIX",
    "VIRTUAL_ENV",
    "OMP_NUM_THREADS",
    "MKL_THREADING_LAYER",
    "KMP_DUPLICATE_LIB_OK",
    "POLARS_MAX_THREADS",
    "POLARS_SKIP_CPU_CHECK",
    "NUMBA_CACHE_DIR",
    "PYVISTA_OFF_SCREEN",
    "HOME",
    "USERPROFILE",
    "OneDrive",
    "OneDriveCommercial",
]

# (label, statement) imported each in a fresh process. Order follows the
# f4enix import chain.
DEPENDENCIES = [
    ("numpy", "import numpy"),
    ("pandas", "import pandas"),
    ("polars", "import polars"),
    ("matplotlib", "import matplotlib"),
    ("matplotlib.pyplot", "import matplotlib.pyplot"),
    ("vtkmodules.vtkCommonCore", "import vtkmodules.vtkCommonCore"),
    ("vtk (full meta-package)", "import vtk"),
    ("pyvista", "import pyvista"),
    ("sklearn", "import sklearn.neighbors"),
    ("tables", "import tables"),
    ("migjorn", "import migjorn"),
    ("pypact", "import pypact"),
    ("actigamma", "import actigamma"),
    ("endf", "import endf"),
    ("docx", "import docx"),
    ("openpyxl", "import openpyxl"),
    ("metadata.version", "import importlib.metadata as m; m.version('f4enix')"),
]

CHILD_TEMPLATE = """
import faulthandler, sys, time
faulthandler.dump_traceback_later({dump}, repeat=True, file=sys.stderr)
_t0 = time.perf_counter()
{stmt}
faulthandler.cancel_dump_traceback_later()
print("__ELAPSED__ %.3f" % (time.perf_counter() - _t0), flush=True)
{after}
"""

THREADPOOL_SNIPPET = """
import json
try:
    import threadpoolctl
    print("__TPINFO__ " + json.dumps(threadpoolctl.threadpool_info()))
except Exception as e:
    print("__TPINFO__ " + json.dumps({"error": repr(e)}))
"""

IMPORTTIME_RE = re.compile(r"^import time:\s*(\d+)\s*\|\s*(\d+)\s*\|(\s*)(\S+)")


class Report:
    def __init__(self, path: str | None) -> None:
        self._fh = open(path, "w", encoding="utf-8") if path else None

    def __call__(self, text: str = "") -> None:
        print(text, flush=True)
        if self._fh:
            self._fh.write(text + "\n")
            self._fh.flush()

    def section(self, title: str) -> None:
        self()
        self("=" * 78)
        self(title)
        self("=" * 78)

    def close(self) -> None:
        if self._fh:
            self._fh.close()


@dataclass
class RunResult:
    name: str
    wall: float
    elapsed: float | None  # time measured inside the child, None if no result
    returncode: int | None
    timed_out: bool
    stdout: str
    stderr: str
    note: str = ""
    importtime: list[tuple[int, int, int, str]] = field(default_factory=list)

    @property
    def status(self) -> str:
        if self.timed_out:
            return "TIMEOUT"
        if self.returncode != 0:
            return f"ERROR({self.returncode})"
        return "ok"

    def faulthandler_dumps(self) -> str:
        lines = [
            ln for ln in self.stderr.splitlines() if not ln.startswith("import time:")
        ]
        return "\n".join(lines).strip()


def run_child(
    python: str,
    name: str,
    stmt: str,
    cwd: str,
    timeout: float,
    dump: int,
    env_set: dict[str, str] | None = None,
    env_unset: list[str] | None = None,
    flags: list[str] | None = None,
    importtime: bool = False,
    after: str = "",
    note: str = "",
) -> RunResult:
    env = dict(os.environ)
    for key in env_unset or []:
        env.pop(key, None)
    env.update(env_set or {})
    cmd = [python] + (flags or [])
    if importtime:
        cmd += ["-X", "importtime"]
    code = CHILD_TEMPLATE.format(dump=dump, stmt=stmt, after=after)
    cmd += ["-c", code]

    t0 = time.perf_counter()
    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    timed_out = False
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        proc.kill()
        out, err = proc.communicate()
    wall = time.perf_counter() - t0

    elapsed = None
    m = re.search(r"__ELAPSED__ ([\d.]+)", out)
    if m:
        elapsed = float(m.group(1))

    res = RunResult(
        name=name,
        wall=wall,
        elapsed=elapsed,
        returncode=None if timed_out else proc.returncode,
        timed_out=timed_out,
        stdout=out,
        stderr=err,
        note=note,
    )
    if importtime:
        for ln in err.splitlines():
            mm = IMPORTTIME_RE.match(ln)
            if mm:
                res.importtime.append(
                    (
                        int(mm.group(1)),
                        int(mm.group(2)),
                        len(mm.group(3)),
                        mm.group(4),
                    )
                )
    return res


# --------------------------------------------------------------------------
# Filesystem detection
# --------------------------------------------------------------------------
def _linux_mounts() -> list[tuple[str, str, str]]:
    mounts = []
    try:
        with open("/proc/mounts", encoding="utf-8") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 3:
                    mnt = parts[1].replace("\\040", " ")
                    mounts.append((mnt, parts[2], parts[0]))
    except OSError:
        pass
    return mounts


NETWORK_FS = {"nfs", "nfs4", "cifs", "smb3", "smbfs", "fuse.sshfs", "gpfs", "lustre",
              "beegfs", "ceph", "fuse.gcsfuse", "9p", "drvfs", "fuse.rclone"}


def fs_info(path: str) -> str:
    """Return a short description of the filesystem holding ``path``."""
    if not path:
        path = os.getcwd()
    try:
        real = os.path.realpath(path)
    except OSError:
        real = path
    if not os.path.exists(real):
        return "missing"

    if IS_LINUX:
        best = ("", "?", "?")
        for mnt, fstype, src in _linux_mounts():
            if (real == mnt or real.startswith(mnt.rstrip("/") + "/")) and len(
                mnt
            ) >= len(best[0]):
                best = (mnt, fstype, src)
        tag = " NETWORK" if best[1] in NETWORK_FS or best[1].startswith("fuse") else ""
        return f"{best[1]} ({best[2]} on {best[0]}){tag}"

    if IS_WIN:
        return _win_fs_info(real)

    return "unknown (unsupported platform)"


def _win_fs_info(path: str) -> str:
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.windll.kernel32
    vol = ctypes.create_unicode_buffer(260)
    if not kernel32.GetVolumePathNameW(path, vol, 260):
        return "unknown"
    root = vol.value
    drive_types = {0: "unknown", 1: "no-root", 2: "removable", 3: "local",
                   4: "NETWORK", 5: "cdrom", 6: "ramdisk"}
    dtype = drive_types.get(kernel32.GetDriveTypeW(root), "?")
    fsname = ctypes.create_unicode_buffer(64)
    kernel32.GetVolumeInformationW(
        root, None, 0, None, None, None, fsname, 64
    )
    tags = []
    if path.startswith("\\\\"):
        tags.append("UNC")
    for var in ("OneDrive", "OneDriveCommercial", "OneDriveConsumer"):
        od = os.environ.get(var)
        if od and os.path.normcase(path).startswith(os.path.normcase(od)):
            tags.append("ONEDRIVE")
            break
    # Files-On-Demand placeholder: FILE_ATTRIBUTE_RECALL_ON_DATA_ACCESS
    attrs = kernel32.GetFileAttributesW(wintypes.LPCWSTR(path))
    if attrs != 0xFFFFFFFF and attrs & 0x00400000:
        tags.append("CLOUD-PLACEHOLDER")
    extra = (" " + " ".join(tags)) if tags else ""
    return f"{fsname.value or '?'} {dtype} ({root}){extra}"


def windows_defender() -> list[str]:
    if not IS_WIN or not shutil.which("powershell"):
        return []
    out = []
    cmds = {
        "RealTimeProtectionEnabled": "(Get-MpComputerStatus).RealTimeProtectionEnabled",
        "ExclusionPath": "(Get-MpPreference).ExclusionPath -join '; '",
        "ExclusionProcess": "(Get-MpPreference).ExclusionProcess -join '; '",
    }
    for label, ps in cmds.items():
        try:
            r = subprocess.run(
                ["powershell", "-NoProfile", "-NonInteractive", "-Command", ps],
                capture_output=True, text=True, timeout=30,
            )
            val = (r.stdout.strip() or r.stderr.strip().splitlines()[0:1] or ["?"])
            out.append(f"{label}: {val if isinstance(val, str) else val[0]}")
        except Exception as exc:  # noqa: BLE001
            out.append(f"{label}: n/a ({exc!r})")
    return out


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def fmt_time(res: RunResult) -> str:
    if res.timed_out:
        return f"> {res.wall:7.1f} s"
    t = res.elapsed if res.elapsed is not None else res.wall
    return f"{t:9.2f} s"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--python", default=sys.executable,
                        help="interpreter to test (default: this one)")
    parser.add_argument("--cwd", default=os.getcwd(),
                        help="working directory for the baseline runs "
                             "(e.g. the notebook folder)")
    parser.add_argument("--timeout", type=float, default=180.0,
                        help="timeout per subprocess in seconds (default 180)")
    parser.add_argument("--dump-interval", type=int, default=20,
                        help="seconds between faulthandler stack dumps (default 20)")
    parser.add_argument("--skip-deps", action="store_true",
                        help="skip per-dependency timings")
    parser.add_argument("--with-dead-display", action="store_true",
                        help="also run with DISPLAY pointing to an unreachable "
                             "X server (Linux; may take up to --timeout)")
    parser.add_argument("--output", help="also write the report to this file")
    args = parser.parse_args()

    rep = Report(args.output)
    py, cwd, to, dump = args.python, os.path.abspath(args.cwd), args.timeout, args.dump_interval
    local_tmp = tempfile.mkdtemp(prefix="f4enix_diag_")
    results: list[RunResult] = []

    # ---------------------------------------------------------------- env
    rep.section("1. ENVIRONMENT")
    rep(f"date            : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    rep(f"platform        : {platform.platform()}  ({platform.machine()})")
    rep(f"cpu count       : {os.cpu_count()}")
    rep(f"tested python   : {py}")
    info = run_child(
        py, "info",
        "import sys, json, platform\n"
        "print('__INFO__ ' + json.dumps({'version': sys.version, 'prefix': sys.prefix,"
        " 'path': sys.path}))",
        cwd=cwd, timeout=60, dump=dump,
    )
    child = {}
    m = re.search(r"__INFO__ (.*)", info.stdout)
    if m:
        child = json.loads(m.group(1))
        rep(f"python version  : {child['version'].splitlines()[0]}")
        rep(f"sys.prefix      : {child['prefix']}  -> {fs_info(child['prefix'])}")
    else:
        rep(f"could not query interpreter: {info.stderr.strip()[:500]}")
    rep(f"baseline cwd    : {cwd}  -> {fs_info(cwd)}")
    rep()
    rep("environment variables:")
    for var in ENV_VARS:
        if var in os.environ:
            rep(f"  {var} = {os.environ[var]}")
    rep()
    rep("sys.path of the tested interpreter when run from baseline cwd "
        "('' = cwd, searched FIRST):")
    for entry in child.get("path", []):
        rep(f"  {entry or repr('')!s:70s} {fs_info(entry or cwd)}")

    # versions without importing (metadata only)
    ver = run_child(
        py, "versions",
        "import importlib.metadata as m, json\n"
        "pk = ['f4enix','migjorn','numjuggler','matplotlib','vtk','pyvista','numpy',"
        "'pandas','polars','scikit-learn','tables','numba','pypact','actigamma','endf',"
        "'threadpoolctl','ipykernel']\n"
        "d = {}\n"
        "for p in pk:\n"
        "    try: d[p] = m.version(p)\n"
        "    except Exception: d[p] = '-'\n"
        "print('__VER__ ' + json.dumps(d))",
        cwd=local_tmp, timeout=120, dump=dump,
    )
    m = re.search(r"__VER__ (.*)", ver.stdout)
    if m:
        rep()
        rep("package versions: " + ", ".join(f"{k}={v}" for k, v in json.loads(m.group(1)).items()))

    mpl = run_child(
        py, "mpl-info",
        "import matplotlib, os, json, glob\n"
        "c = matplotlib.get_cachedir()\n"
        "print('__MPL__ ' + json.dumps({'cachedir': c, 'configdir': matplotlib.get_configdir(),"
        " 'writable': os.access(c, os.W_OK),"
        " 'fontlists': sorted(os.path.basename(f) for f in glob.glob(os.path.join(c, 'fontlist-*'))),"
        " 'locks': glob.glob(os.path.join(c, '*.matplotlib-lock'))}))",
        cwd=local_tmp, timeout=to, dump=dump,
    )
    m = re.search(r"__MPL__ (.*)", mpl.stdout)
    rep()
    if m:
        d = json.loads(m.group(1))
        rep(f"matplotlib cache: {d['cachedir']} (writable={d['writable']}) -> {fs_info(d['cachedir'])}")
        rep(f"  font lists    : {d['fontlists'] or 'NONE (will be rebuilt on next pyplot import)'}")
        if d["locks"]:
            rep(f"  STALE LOCKS   : {d['locks']}")
    else:
        rep(f"matplotlib info failed ({mpl.status}): {mpl.stderr.strip()[-500:]}")

    if IS_WIN:
        rep()
        rep("Windows Defender (best effort, exclusions may require admin):")
        for line in windows_defender():
            rep(f"  {line}")

    # ---------------------------------------------------------- scenarios
    rep.section("2. A/B SCENARIOS: import f4enix")
    fresh_mpl = tempfile.mkdtemp(prefix="f4enix_diag_mpl_")
    scenarios = [
        dict(name="baseline (cold?)", cwd=cwd, importtime=True,
             note="first run: includes cold disk cache / antivirus scan"),
        dict(name="baseline (repeat)", cwd=cwd,
             note="big gap vs first run -> cold cache / antivirus"),
        dict(name="cwd off sys.path (-P, local tmp cwd)", cwd=local_tmp, flags=["-P"],
             note="much faster -> cwd/network share is the problem"),
        dict(name="MPLBACKEND=Agg", cwd=cwd, env_set={"MPLBACKEND": "Agg"},
             note="faster -> matplotlib backend detection"),
        dict(name="fresh matplotlib cache", cwd=cwd, env_set={"MPLCONFIGDIR": fresh_mpl},
             note="forces font-cache rebuild (cost a new env pays once)"),
    ]
    if not IS_WIN:
        scenarios.append(dict(name="DISPLAY unset", cwd=cwd,
                              env_unset=["DISPLAY", "WAYLAND_DISPLAY"],
                              note="faster -> X display probing"))
        if args.with_dead_display:
            scenarios.append(dict(name="DISPLAY=10.255.255.1:0 (dead)", cwd=cwd,
                                  env_set={"DISPLAY": "10.255.255.1:0"},
                                  env_unset=["MPLBACKEND"],
                                  note="reproduces a stale X11 forwarding"))
    for sc in scenarios:
        rep(f"running: {sc['name']} ...")
        res = run_child(py, stmt="import f4enix", timeout=to, dump=dump,
                        after=THREADPOOL_SNIPPET if sc["name"].startswith("baseline (repeat") else "",
                        **sc)
        results.append(res)
        rep(f"  -> {res.status:10s} {fmt_time(res)}")

    rep()
    rep(f"{'scenario':42s} {'status':10s} {'time':>11s}  hint")
    for res in results:
        rep(f"{res.name:42s} {res.status:10s} {fmt_time(res)}  {res.note}")

    # ------------------------------------------------------------- deps
    dep_results: list[RunResult] = []
    if not args.skip_deps:
        rep.section("3. PER-DEPENDENCY IMPORT TIME (each in a fresh process, cwd = baseline)")
        for label, stmt in DEPENDENCIES:
            res = run_child(py, label, stmt, cwd=cwd, timeout=to, dump=dump)
            dep_results.append(res)
            err = ""
            if res.status.startswith("ERROR"):
                last = [ln for ln in res.stderr.strip().splitlines() if ln][-1:]
                err = last[0][:80] if last else ""
            rep(f"  {label:32s} {res.status:10s} {fmt_time(res)}  {err}")

    # --------------------------------------------------------- openmp
    rep.section("4. OPENMP / BLAS RUNTIMES LOADED AFTER import f4enix")
    tp = next((r for r in results if r.name == "baseline (repeat)"), None)
    m = re.search(r"__TPINFO__ (.*)", tp.stdout) if tp else None
    if m:
        info_tp = json.loads(m.group(1))
        if isinstance(info_tp, dict):
            rep(f"  threadpoolctl unavailable: {info_tp.get('error')}")
        else:
            omp = set()
            for lib in info_tp:
                rep(f"  {lib.get('user_api'):8s} {lib.get('internal_api'):10s} "
                    f"threads={lib.get('num_threads')}  {lib.get('filepath')}")
                if lib.get("user_api") == "openmp":
                    omp.add(lib.get("internal_api"))
            if len(omp) > 1:
                rep(f"  WARNING: several OpenMP runtimes loaded {sorted(omp)}: "
                    "known source of hangs/deadlocks")
    else:
        rep("  not available (baseline run failed?)")

    # --------------------------------------------------------- details
    rep.section("5. DETAILS")
    base = results[0]
    if base.importtime:
        rep("slowest modules of the first baseline run (-X importtime, cumulative):")
        top = sorted(base.importtime, key=lambda x: x[1], reverse=True)[:15]
        for self_us, cum_us, _indent, mod in top:
            rep(f"  {cum_us / 1e6:8.2f} s cumulative  {self_us / 1e6:8.2f} s self   {mod}")
        rep("  (self time = time spent in that module itself, excluding its imports)")
    for res in results + dep_results:
        if res.timed_out or res.returncode != 0:
            rep()
            rep(f"--- {res.name}: {res.status} ---")
            dumps = res.faulthandler_dumps()
            rep(dumps[-6000:] if dumps else "(no stderr)")

    rep.section("DONE")
    rep("Please share this whole output. See the docstring of this script for "
        "how to get a stack trace from inside a hanging Jupyter kernel.")
    rep.close()
    shutil.rmtree(fresh_mpl, ignore_errors=True)
    shutil.rmtree(local_tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
