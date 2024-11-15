"""
Auxiliary functions used by different modules

"""

from __future__ import annotations

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

import os
import re

from numjuggler.parser import Card

from f4enix.constants import PAT_COMMENT_TEXT, PAT_DOLLAR_COMMENT

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
    comments = ""
    for line in card.lines:
        # check if either dollar or c comments are matched and store them
        c_comm = PAT_COMMENT_TEXT.match(line)
        d_comm = PAT_DOLLAR_COMMENT.search(line)

        for match in [c_comm, d_comm]:
            if match is not None:
                comments = comments + match.group()

    return comments


def _detect_decoding_errors_line(l, _s=_surrogates.finditer):  # pragma: no cover
    """Return decoding errors in a line of text

    Works with text lines decoded with the surrogateescape
    error handler.

    Returns a list of (pos, byte) tuples

    """
    # DC80 - DCFF encode bad bytes 80-FF
    return [(m.start(), bytes([ord(m.group()) - 0xDC00])) for m in _s(l)]


def debug_file_unicode(file: os.PathLike | str) -> str:  # pragma: no cover
    """given a file prints the unicode error that were found and at what
    line

    Parameters
    ----------
    file : os.PathLike | str
        file to be debugged

    Returns
    -------
    str
        bug found
    """
    txt = ""
    with open(file, "r", errors="surrogateescape") as f:
        for i, line in enumerate(f, 1):
            errors = _detect_decoding_errors_line(line)
            if errors:
                txt += f"Found errors on line {i}:"
                for col, b in errors:
                    txt += f" {col + 1:2d}: {b[0]:02x}"
                txt += "\n"
    return txt
