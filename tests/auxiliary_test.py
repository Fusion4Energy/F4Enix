from importlib.resources import files, as_file
from numjuggler.parser import Card
import pytest

from f4enix.input.auxiliary import get_comments


@pytest.mark.parametrize(
    ["lines", "comment"],
    [[["C dasd\n", "M10 8016 1 $inline comment\n"], "C dasd\n$inline comment\n"]],
)
def test_get_comments(lines, comment):
    card = Card(lines, 5, 10)
    card.get_values()
    comments = get_comments(card)
    assert comments == comment
