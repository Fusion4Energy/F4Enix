import migjorn
from collections.abc import MutableMapping
import logging


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
        return (str(c.id) for c in self._model.cells())

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
        return (self._key(s) for s in self._model.surfaces())

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
        try:
            tid = int(key.upper().lstrip("*").removeprefix("TR"))
        except ValueError:
            raise KeyError(key)
        text = value.text if isinstance(value, migjorn.Transform) else value
        self._model.remove_transform(tid)
        self._model.add_transform(text)

    def __delitem__(self, key: str) -> None:
        try:
            tid = int(key.upper().lstrip("*").removeprefix("TR"))
        except ValueError:
            raise KeyError(key)
        self._model.remove_transform(tid)

    def __iter__(self):
        return (f"TR{t.id}" for t in self._model.transforms())

    def __len__(self) -> int:
        return self._model.num_transforms


class _OtherDataProxy(MutableMapping):
    """Dict-like proxy for other_data; delegates directly to the underlying dict."""

    def _match(self, key: str, card: migjorn.DataCard) -> bool:
        """Check if the key matches the card name, handles particles (e.g. F6:N,P)"""
        name = card.name
        if name is None:
            return False

        particles = None
        key = key.lower()

        if ":" in key:
            key, particles = key.split(":")

        if key == name.lower():
            if particles is None:
                return True
            else:
                if card.particle.lower() == particles:
                    return True
        return False

    @staticmethod
    def _key(card: migjorn.DataCard) -> str:
        return str(card.name) + (f":{card.particle}" if card.particle else "")

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    def __getitem__(self, key: str) -> migjorn.DataCard:
        # check if particles are specified in the key (e.g. F6:N,P)

        for card in self._model.data_cards():
            if self._match(key, card):
                return card

        raise KeyError(f"Card {key} not found in data cards")

    def __setitem__(self, key: str, value: str) -> None:
        try:
            del self[key]  # remove existing card if it exists
        except KeyError:
            pass  # it is ok, remove it only if found
        self._model.add_data_card(value)

    def __delitem__(self, key: str) -> None:
        for card in self._model.data_cards():
            if self._match(key, card):
                card.remove()
                return
        raise KeyError(f"Card {key} not found in data cards")

    def __iter__(self):
        # exclude materials and transforms
        # TODO this will need to be changed from migjorn
        keys = []
        none_cards = 0
        for card in self._model.data_cards():
            if card.name is None:
                name = f"NONE{none_cards}"
                none_cards += 1
            else:
                name = card.name
            if not (name.lower().startswith("m") or name.lower().startswith("tr")):
                keys.append(self._key(card))
        return iter(keys)

    def __len__(self) -> int:
        return len(self._model.data_cards())
