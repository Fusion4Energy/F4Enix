from importlib.resources import as_file, files

import pytest

from f4enix.core.auxiliary import get_comments


@pytest.mark.parametrize(
    ["text", "comment"],
    [["C dasd\nM10 8016 1 $inline comment\n", "C dasd\n$inline comment\n"]],
)
def test_get_comments(text, comment):
    comments = get_comments(text)
    assert comments == comment
