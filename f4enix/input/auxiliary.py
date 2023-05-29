"""
Auxiliary functions used by different modules

"""
from numjuggler.parser import Card
from f4enix.constants import PAT_COMMENT, PAT_DOLLAR_COMMENT, PAT_COMMENT_TEXT
import re
import os

_surrogates = re.compile(r"[\uDC80-\uDCFF]")


def get_comments(card: Card) -> str:
    """Analyze the lines of a card and get the comments

    Parameters
    ----------
    card : Card
        card to be analyzed

    Returns
    -------
    str
        comments as a single string
    """
    comments = ''
    for line in card.lines:
        # check if either dollar or c comments are matched and store them
        c_comm = PAT_COMMENT_TEXT.match(line)
        d_comm = PAT_DOLLAR_COMMENT.search(line)

        for match in [c_comm, d_comm]:
            if match is not None:
                comments = comments+match.group()

    return comments


def _detect_decoding_errors_line(l, _s=_surrogates.finditer):  # pragma: no cover
    """Return decoding errors in a line of text

    Works with text lines decoded with the surrogateescape
    error handler.

    Returns a list of (pos, byte) tuples

    """
    # DC80 - DCFF encode bad bytes 80-FF
    return [(m.start(), bytes([ord(m.group()) - 0xDC00]))
            for m in _s(l)]


def debug_file_unicode(file: os.PathLike) -> str:  # pragma: no cover
    """given a file prints the unicode error that were found and at what
    line

    Parameters
    ----------
    file : os.PathLike
        file to be debugged

    Returns
    -------
    str
        bug found
    """
    txt = ''
    with open(file, 'r', errors="surrogateescape") as f:
        for i, line in enumerate(f, 1):
            errors = _detect_decoding_errors_line(line)
            if errors:
                txt += (f"Found errors on line {i}:")
                for (col, b) in errors:
                    txt += (f" {col + 1:2d}: {b[0]:02x}")
                txt += '\n'
    return txt
