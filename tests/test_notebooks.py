import os
import sys
from pathlib import Path

import nbformat
import pytest
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert.preprocessors.execute import CellExecutionError


def _notebook_run(path):
    """
    Execute a notebook via nbconvert and collect output.
    :returns (parsed nb object, execution errors)
    """
    kernel_name = "python%d" % sys.version_info[0]
    errors = []

    with open(path, encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)
        nb.metadata.get("kernelspec", {})["name"] = kernel_name
        ep = ExecutePreprocessor(
            kernel_name=kernel_name, timeout=300
        )  # , allow_errors=True

        try:
            ep.preprocess(nb, {"metadata": {"path": path.parent}})
            print(f"running notebook from this path {path.parent}")

        except CellExecutionError as e:
            if "SKIP" in e.traceback:
                print(str(e.traceback).split("\n")[-2])
            else:
                raise e

    return nb, errors


WINDOWS_ACCESS_ERROR = ["meshtal.ipynb", "tutorial.ipynb"]


@pytest.mark.parametrize(
    "filename", Path(os.path.join("docs", "source")).rglob("*.ipynb")
)
def test_task_1(filename):
    if os.path.basename(filename) in WINDOWS_ACCESS_ERROR:
        return

    print(f"Attempting to run {filename}")
    _, errors = _notebook_run(filename)
    assert errors == []


@pytest.mark.skipif(sys.platform == "win32", reason="Windows access error")
def test_plotting_jupyters():
    _, errors = _notebook_run(Path("docs/source/tutorial/tutorial.ipynb"))
    assert errors == []

    _, errors = _notebook_run(
        Path("docs/source/examples/output/jupyters/meshtal.ipynb")
    )
    assert errors == []
