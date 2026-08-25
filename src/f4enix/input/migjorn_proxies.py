import migjorn
from collections.abc import MutableMapping
import re
from f4enix.core.errors import InvalidCardError

PAT_NOT_OTHER = re.compile(r"(m|mx|tr)\d+", re.IGNORECASE)


class _CellsProxy(MutableMapping[str, migjorn.Cell]):
    """Dict-like proxy for model cells; __setitem__ replaces a cell in the model.

    Attributes
    ----------
    _model : migjorn.Model
        the underlying migjorn model whose cells are proxied.
    """

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
        cell = self._model.add_cell(text)
        if not cell.well_formed:
            raise InvalidCardError(cell)

    def __delitem__(self, key: str) -> None:
        self._model.remove_cell(int(key))

    def __iter__(self):
        return (str(c.id) for c in self._model.cells())

    def __len__(self) -> int:
        return self._model.num_cells


class _SurfsProxy(MutableMapping[str, migjorn.Surface]):
    """Dict-like proxy for model surfaces; __setitem__ replaces a surface in the model.

    Attributes
    ----------
    _model : migjorn.Model
        the underlying migjorn model whose surfaces are proxied.
    """

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
        sur = self._model.add_surface(text)
        if not sur.well_formed:
            raise InvalidCardError(sur)

    def __delitem__(self, key: str) -> None:
        self._model.remove_surface(int(key.lstrip("*")))

    def __iter__(self):
        return (self._key(s) for s in self._model.surfaces())

    def __len__(self) -> int:
        return self._model.num_surfaces


class _TransformsProxy(MutableMapping[str, migjorn.Transform]):
    """Dict-like proxy for model transforms keyed by 'TRn'.

    Deletion is supported via ``del inp.transformations['TR5']``.
    Setting is not supported (migjorn has no add_transform).

    Attributes
    ----------
    _model : migjorn.Model
        the underlying migjorn model whose transforms are proxied.
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


class _OtherDataProxy(MutableMapping[str, migjorn.DataCard]):
    """Dict-like proxy for other_data; delegates directly to the underlying dict.

    Attributes
    ----------
    _model : migjorn.Model
        the underlying migjorn model whose data cards are proxied.
    """

    def _match(
        self, key: str, card: migjorn.DataCard, none_index: int | None
    ) -> bool:
        """Check if the key matches the card, handles particles (e.g. F6:N,P)

        Cards with no parseable name (e.g. bare $-comment lines or informal
        tables living in the data block) are matched against their ordinal
        placeholder key instead (e.g. 'NONE0').
        """
        name = card.name
        if name is None:
            return none_index is not None and key == f"NONE{none_index}"

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
    def _key(card: migjorn.DataCard, none_index: int | None = None) -> str:
        if card.name is None:
            return f"NONE{none_index}"
        return card.name + (f":{card.particle}" if card.particle else "")

    def __init__(self, model: migjorn.Model) -> None:
        self._model = model

    def _iter_cards(self):
        """Yield (card, none_index) for every data card, numbering the nameless
        ones in iteration order so they get a stable, distinct placeholder key."""
        none_idx = 0
        for card in self._model.data_cards():
            if card.name is None:
                yield card, none_idx
                none_idx += 1
            else:
                yield card, None

    def __getitem__(self, key: str) -> migjorn.DataCard:
        # check if particles are specified in the key (e.g. F6:N,P)

        for card, none_idx in self._iter_cards():
            if self._match(key, card, none_idx):
                return card

        raise KeyError(f"Card {key} not found in data cards")

    def __setitem__(self, key: str, value: str | migjorn.DataCard) -> None:
        try:
            del self[key]  # remove existing card if it exists
        except KeyError:
            pass  # it is ok, remove it only if found
        if isinstance(value, str):
            self._model.add_data_card(value)
        else:
            self._model.add_data_card(value.text)

    def __delitem__(self, key: str) -> None:
        for card, none_idx in self._iter_cards():
            if self._match(key, card, none_idx):
                card.remove()
                return
        raise KeyError(f"Card {key} not found in data cards")

    def __iter__(self):
        # exclude materials and transforms, keep cards with no parseable name
        # (e.g. bare $-comment lines or informal tables in the data block)
        # TODO this will need to be changed from migjorn
        keys = []
        for card, none_idx in self._iter_cards():
            if card.name is None:
                keys.append(self._key(card, none_idx))
            elif not PAT_NOT_OTHER.match(card.name):
                keys.append(self._key(card))
        return iter(keys)

    def __len__(self) -> int:
        return sum(1 for _ in self)
