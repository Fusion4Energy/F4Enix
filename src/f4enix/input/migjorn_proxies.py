import migjorn
from collections.abc import MutableMapping

from f4enix.core.constants import PAT_F_TR_CARD_KEY

# TODO: migjorn has no attribute text for datacards


class _CellsProxy(MutableMapping):
    """Dict-like proxy for model cells; __setitem__ replaces a cell in the model."""

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    def __getitem__(self, key: str) -> migjorn.Cell:
        cell = self._model.cell(int(key))
        if cell is None:
            raise KeyError(key)
        return cell

    def __setitem__(self, key: str, value: "migjorn.Cell | str") -> None:
        """Replace the cell in the model. value may be a Cell handle or raw text."""
        cell_id = int(key)
        text = value.text if isinstance(value, migjorn.Cell) else value
        self._model.remove_cell(cell_id)
        self._model.add_cell(text)

    def __delitem__(self, key: str) -> None:
        self._model.remove_cell(int(key))

    def __iter__(self):
        return (str(c.id) for c in self._model.cells)

    def __len__(self) -> int:
        return self._model.num_cells


class _SurfsProxy(MutableMapping):
    """Dict-like proxy for model surfaces; __setitem__ replaces a surface in the model."""

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    @staticmethod
    def _key(surf: migjorn.Surface) -> str:
        return ("*" if surf.reflective else "") + str(surf.id)

    def __getitem__(self, key: str) -> migjorn.Surface:
        sid = int(key.lstrip("*"))
        surf = self._model.surface(sid)
        if surf is None:
            raise KeyError(key)
        return surf

    def __setitem__(self, key: str, value: "migjorn.Surface | str") -> None:
        """Replace the surface in the model. value may be a Surface handle or raw text."""
        sid = int(key.lstrip("*"))
        text = value.text if isinstance(value, migjorn.Surface) else value
        self._model.remove_surface(sid)
        self._model.add_surface(text)

    def __delitem__(self, key: str) -> None:
        self._model.remove_surface(int(key.lstrip("*")))

    def __iter__(self):
        return (self._key(s) for s in self._model.surfaces)

    def __len__(self) -> int:
        return self._model.num_surfaces


class _TransformsProxy(MutableMapping):
    """Dict-like proxy for model transforms keyed by 'TRn'.

    Deletion is supported via ``del inp.transformations['TR5']``.
    Setting is not supported (migjorn has no add_transform).
    """

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    def __getitem__(self, key: str) -> migjorn.Transform:
        try:
            tid = int(key.upper().lstrip("*").removeprefix("TR"))
        except ValueError:
            raise KeyError(key)
        t = self._model.transform(tid)
        if t is None:
            raise KeyError(key)
        return t

    def __setitem__(self, key: str, value) -> None:
        raise NotImplementedError(
            "migjorn does not expose add_transform; modify the transform handle in-place."
        )

    def __delitem__(self, key: str) -> None:
        try:
            tid = int(key.upper().lstrip("*").removeprefix("TR"))
        except ValueError:
            raise KeyError(key)
        self._model.remove_transform(tid)

    def __iter__(self):
        return (f"TR{t.id}" for t in self._model.transforms)

    def __len__(self) -> int:
        return self._model.num_transforms


class _OtherDataProxy(MutableMapping):
    """Dict-like proxy for other_data; delegates directly to the underlying dict."""

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    def __getitem__(self, key: str) -> migjorn.DataCard:
        # check if particles are specified in the key (e.g. F6:N,P)
        if ":" in key:
            key, particles = key.split(":")
        else:
            particles = None
        for card in self._model.data_cards:
            if key == card.name:
                if particles is None:
                    return card
                else:
                    if card.particle == particles.lower():
                        return card

        raise KeyError(f"Card {key} not found in data cards")

    def __setitem__(self, key: str, value: str) -> None:
        # TODO: there is no setter
        raise NotImplementedError()

    def __delitem__(self, key: str) -> None:
        # TODO: there is no simple way to delete a data card
        raise NotImplementedError()

    def __iter__(self):
        return iter(self._model.data_cards)

    def __len__(self) -> int:
        return len(self._model.data_cards)
