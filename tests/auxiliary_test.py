from importlib.resources import as_file, files

import pytest
from numjuggler.parser import Card

from f4enix.core.auxiliary import get_comments


@pytest.mark.parametrize(
    ["lines", "comment"],
    [[["C dasd\n", "M10 8016 1 $inline comment\n"], "C dasd\n$inline comment\n"]],
)
def test_get_comments(lines, comment):
    card = Card(lines, 5, 10)
    card.get_values()
    comments = get_comments(card)
    assert comments == comment
