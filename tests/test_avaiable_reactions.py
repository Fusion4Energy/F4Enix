from f4enix.core.available_reactions import AVAILABLE_REACTIONS


def test_available_reactions():
    for _, pathways in AVAILABLE_REACTIONS.items():
        # check that there are no duplicates
        reactions = [p.__str__() for p in pathways.pathways]
        assert len(reactions) == len(set(reactions))
