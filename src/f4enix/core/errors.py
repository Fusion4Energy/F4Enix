from migjorn import Cell, Surface, Transform


class InvalidCardError(Exception):
    """Exception raised for invalid MCNP cards."""

    def __init__(self, card: Cell | Surface | Transform | str) -> None:
        self.card = card if isinstance(card, str) else card.text
        self.message = "Invalid MCNP card"
        super().__init__(f"{self.message}: {self.card}")
