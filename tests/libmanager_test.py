from importlib.resources import as_file, files

import pandas as pd
import pytest

import f4enix.resources as pkg_res
import tests.resources.libmanager as lib_res
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import Zaid

resources = files(pkg_res)
lib_resources = files(lib_res)
ACTIVATION_FILE = as_file(lib_resources.joinpath("Activation libs.xlsx"))
XSDIR_FILE = as_file(lib_resources.joinpath("xsdir"))
ISOTOPES_FILE = as_file(resources.joinpath("Isotopes.txt"))


class TestLibManger:
    with (
        XSDIR_FILE as xsdir_file,
        ACTIVATION_FILE as activation_file,
        ISOTOPES_FILE as isotopes_file,
    ):
        lm = LibManager(
            xsdir_file, activationfile=activation_file, isotopes_file=isotopes_file
        )

    def test_reactionfilereading(self):
        assert len(self.lm.reactions["99c"]) == 100
        assert len(self.lm.reactions["98c"]) == 34

    def test_default_lm(self):
        LibManager()
        assert True

    def test_get_reactions1(self):
        """
        Test ability to recover reactions for parent zaid (one)
        """
        parent = "9019"
        reaction = self.lm.get_reactions("99c", parent)[0]
        assert reaction[0] == "16"
        assert reaction[1] == "9018"

    def test_get_reactions2(self):
        """
        Test ability to recover reactions for parent zaid (multiple)
        """
        parent = "11023"
        reactions = self.lm.get_reactions("99c", parent)
        print(reactions)
        reaction1 = reactions[0]
        reaction2 = reactions[1]

        assert reaction1[0] == "16"
        assert reaction1[1] == "11022"
        assert reaction2[0] == "102"
        assert reaction2[1] == "11024"

    def test_get_reactions_default(self):
        lm = LibManager()
        # test for metastable reactions
        parent = "78195"  # Pt-195
        reaction = lm.get_reactions("93c", parent)[0]
        assert reaction[0] == "304"
        assert reaction[1] == "78195900"

    def test_formula_conversion(self):
        """
        Test the abilty to switch between isotopes formulas and zaid number
        """
        tests = ["N15", "Er164", "Kr83"]
        finals = ["N-15", "Er-164", "Kr-83"]
        for test, final in zip(tests, finals):
            conversion = self.lm.get_zaidnum(test)
            name, formula = self.lm.get_zaidname(conversion)

            assert final == formula

    def test_check4zaid(self):
        """
        Correctly checks availability of zaids
        """
        zaid = "1001"
        libs = self.lm.check4zaid(zaid)
        assert len(libs) > 1
        assert len(libs[0]) == 3

        zaid = "1010"
        assert len(self.lm.check4zaid(zaid)) == 0

    def test_convertZaid(self):
        # --- Exception if library is not available ---
        try:
            zaid = "1001"
            lib = "44c"
            self.lm.convertZaid(zaid, lib)
            assert False
        except ValueError:
            assert True

        # --- Natural zaid ---
        # 1 to 1
        zaid = "12000"
        lib = "21c"
        translation = self.lm.convertZaid(zaid, lib)
        assert translation == {zaid: (lib, 1, 1)}
        # expansion
        lib = "31c"
        translation = self.lm.convertZaid(zaid, lib)
        assert len(translation) == 3
        # not available in the requested lib but available in default

        # not available
        try:
            zaid = "84000"
            translation = self.lm.convertZaid(zaid, lib)
            assert False
        except ValueError:
            assert True

        # --- 1 to 1 ---
        zaid = "1001"
        translation = self.lm.convertZaid(zaid, lib)
        assert translation == {zaid: (lib, 1, 1)}

        # --- absent ---
        # Use the natural zaid
        zaid = "12024"
        lib = "21c"
        translation = self.lm.convertZaid(zaid, lib)
        assert translation == {"12000": (lib, 1, 1)}

        # zaid available in default or other library
        zaid = "84210"
        translation = self.lm.convertZaid(zaid, lib)
        assert translation[zaid][0] != lib

        # zaid does not exist or not available in any library
        zaid = "84200"
        try:
            translation = self.lm.convertZaid(zaid, lib)
            assert False
        except ValueError:
            assert True

    def test_get_libzaids(self):
        lib = "44c"
        zaids = self.lm.get_libzaids(lib)
        assert len(zaids) == 0

        lib = "21c"
        zaids = self.lm.get_libzaids(lib)
        assert len(zaids) == 76
        assert zaids[0] == "1001"

    def test_get_zaidname(self):
        zaid = "1001"
        name, formula = self.lm.get_zaidname(zaid)
        assert name == "hydrogen"
        assert formula == "H-1"

        zaid = "1000"
        name, formula = self.lm.get_zaidname(zaid)
        assert name == "hydrogen"
        assert formula == "H"

    def test_get_zaidnum(self):
        zaid = "92235"
        try:
            zaidnum = self.lm.get_zaidnum(zaid)
            assert False
        except ValueError:
            assert True

        zaid = "U235"
        zaidnum = self.lm.get_zaidnum(zaid)
        assert zaidnum == "92235"

        # check natural zaids
        assert "1000" == self.lm.get_zaidnum("H")

        with pytest.raises(KeyError):
            self.lm.get_zaidnum("Cj")

    def test_get_zaid_mass(self):
        # Normal zaid
        zaid = "99235.31c  -1"
        zaid = Zaid.from_string(zaid)
        mass = self.lm.get_zaid_mass(zaid)
        assert mass == 252.08298
        # Natural zaid
        zaid = "8000.21c  1"
        zaid = Zaid.from_string(zaid)
        mass = self.lm.get_zaid_mass(zaid)
        assert mass == 15.99937442590581

    def test_default_file_generation(self):
        # both default files
        lm = LibManager()

        # these need to be reinitialized because they have already been
        # "consumed"
        resources2 = files(pkg_res)
        XSDIR_FILE2 = as_file(resources2.joinpath("xsdir.txt"))
        ISOTOPES_FILE2 = as_file(resources2.joinpath("Isotopes.txt"))

        # only xsdir default
        with ISOTOPES_FILE2 as isotopes_file:
            lm = LibManager(isotopes_file=isotopes_file)
        # only isotopes default
        with XSDIR_FILE2 as xsdir_file:
            lm = LibManager(xsdir_path=xsdir_file)


class TestMultiCodeLibManger:
    @pytest.fixture
    def lm(self):
        resources2 = files(pkg_res)
        XSDIR_FILE2 = as_file(resources2.joinpath("xsdir.txt"))
        with XSDIR_FILE2 as xsdir_file:
            df_rows = [
                ["99c", "sda", "", xsdir_file],
                ["98c", "acsdc", "", xsdir_file],
                ["21c", "adsadsa", "", xsdir_file],
                ["31c", "adsadas", "", xsdir_file],
                ["00c", "sdas", "yes", xsdir_file],
                ["71c", "sdasxcx", "", xsdir_file],
            ]
        df_lib = pd.DataFrame(df_rows)
        df_lib.columns = ["Suffix", "Name", "Default", "MCNP"]
        resources2 = files(pkg_res)
        lib_resources2 = files(lib_res)
        ACT_FILE2 = as_file(lib_resources2.joinpath("Activation libs.xlsx"))
        ISOTOPES_FILE2 = as_file(resources2.joinpath("Isotopes.txt"))
        with ISOTOPES_FILE2 as isotopes_file, ACT_FILE2 as act_file:
            lm = LibManager(
                df_lib, isotopes_file=isotopes_file, activationfile=act_file
            )

        return lm

    def test_available_libs(self, lm: LibManager):
        """
        Test the ability to recover the available libraries
        """
        libs = lm.libraries["mcnp"]
        assert "99c" in libs
        assert "98c" in libs
        assert "21c" in libs
        assert "31c" in libs
        assert "00c" in libs
        assert "71c" in libs
