"""
Tests to prevent regressions in import time.

The root cause of the BUG "Import time is sometimes 10+ minutes":
  - ``MCNPinput.py`` imports ``matplotlib.pyplot`` at module level.
  - This triggers ``matplotlib.font_manager`` which scans **all system fonts**
    when its JSON cache is cold (first run after install or after a matplotlib
    upgrade).
  - On machines with many fonts (Windows with Office/Adobe fonts, Linux HPC
    clusters with hundreds of installed fonts) this scan can take 10+ minutes.

The invariant we enforce:
  importing ``f4enix`` must NOT trigger ``matplotlib.font_manager``
  (which would force a font-cache rebuild on first import).
"""

import subprocess
import sys


def _run_import_check(code: str) -> subprocess.CompletedProcess:
    """Run *code* in a fresh Python interpreter and return the result."""
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )


def test_import_does_not_trigger_matplotlib_font_manager():
    """Importing f4enix must not load ``matplotlib.font_manager``.

    ``matplotlib.font_manager`` initialises a ``FontManager`` instance at
    module level (``fontManager = _load_fontmanager()``).  When the font
    cache is missing that call scans every font file on the system, which is
    extremely slow on machines with many fonts.

    After the fix, none of the modules imported by ``import f4enix`` should
    pull in ``matplotlib.font_manager`` as a side-effect.
    """
    code = """
import sys

# snapshot sys.modules before the import
before = set(sys.modules)

import f4enix

after = set(sys.modules)
new_modules = after - before

# matplotlib.font_manager is the slow module – it must not be imported eagerly
if "matplotlib.font_manager" in new_modules:
    print("FAIL: matplotlib.font_manager was loaded during 'import f4enix'")
    sys.exit(1)
else:
    print("PASS")
    sys.exit(0)
"""
    result = _run_import_check(code)
    assert result.returncode == 0, (
        "matplotlib.font_manager was eagerly loaded during 'import f4enix'.\n"
        "This causes slow import times on machines with many installed fonts "
        "(the font cache must be rebuilt from scratch).\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr[:400]}"
    )
