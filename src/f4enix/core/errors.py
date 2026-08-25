from migjorn import Cell, Surface


class InvalidCardError(Exception):
    """Exception raised for invalid MCNP cards."""

    def __init__(self, card: Cell | Surface) -> None:
        self.card = card.text
        self.message = "Invalid MCNP card"
        super().__init__(f"{self.message}: {self.card}")
