from numjuggler.parser import Card
from f4eparser.constants import PAT_COMMENT, PAT_DOLLAR_COMMENT


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
        c_comm = PAT_COMMENT.match(line)
        d_comm = PAT_DOLLAR_COMMENT.search(line)

        for match in [c_comm, d_comm]:
            if match is not None:
                comments = comments+match

    return comments
