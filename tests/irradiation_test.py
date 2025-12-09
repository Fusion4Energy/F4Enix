from f4enix.input.irradiation import Nuclide


class TestNuclide:
    def test_from_formula_basic(self):
        n = Nuclide.from_formula("Li3")
        assert n.zaid == 3003
        assert not n.metastable
        assert not n.IRS_active
        assert n.lib is None

    def test_from_formula_complex(self):
        n = Nuclide.from_formula("irsLi3m.99c")
        assert n.zaid == 3003
        assert n.metastable
        assert n.IRS_active
        assert n.lib == "99c"

    def test_from_int_string_basic(self):
        n = Nuclide.from_int_string("3003")
        assert n.zaid == 3003
        assert not n.metastable
        assert not n.IRS_active
        assert n.lib is None

    def test_from_int_string_complex(self):
        n = Nuclide.from_int_string("9993003900.99c")
        assert n.zaid == 3003
        assert n.metastable
        assert n.IRS_active
        assert n.lib == "99c"

    def test_write_to_formula(self):
        n = Nuclide(3003, metastable=True, IRS_active=True, lib="99c")
        formula = n.write_to_formula()
        assert formula == "irsLi3m.99c"

    def test_write_to_int_string(self):
        n = Nuclide(3003, metastable=True, IRS_active=True, lib="99c")
        int_string = n.write_to_int_string()
        assert int_string == "9993003900.99c"
