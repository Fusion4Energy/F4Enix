import os
import re
from copy import deepcopy
from importlib.resources import as_file, files

import pytest

import f4enix.resources as pkg_res
import tests.resources.materials as mat_res
from f4enix.input.libmanager import LibManager
from f4enix.input.materials import Element, MatCardsList, Material, SubMaterial, Zaid
from f4enix.input.MCNPinput import Input

resources = files(mat_res)
resources_pkg = files(pkg_res)

XSDIR = as_file(resources.joinpath("xsdir_mcnp6.2"))
ISOTOPES_FILE = as_file(resources_pkg.joinpath("Isotopes.txt"))
# Other
with XSDIR as xsdir_file, ISOTOPES_FILE as isotopes_file:
    LIBMAN = LibManager(xsdir_file, isotopes_file=isotopes_file, defaultlib="81c")

# Files
with as_file(resources.joinpath("mat_test.i")) as inp:
    inp_matcard1 = MatCardsList.from_input(inp)
    inp_matcard2 = Input.from_input(inp).materials

with as_file(resources.joinpath("mat_test2.i")) as inp2:
    inp2_matcard1 = MatCardsList.from_input(inp2)
    inp2_matcard2 = Input.from_input(inp2).materials

with as_file(resources.joinpath("test.i")) as inp3:
    inp3_matcard1 = MatCardsList.from_input(inp3)
    inp3_matcard2 = Input.from_input(inp3).materials

with as_file(resources.joinpath("activation.i")) as inp:
    activation_matcard1 = MatCardsList.from_input(inp)
    activation_matcard2 = Input.from_input(inp).materials


class TestZaid:
    tests = [
        {"str": "1001.31c   -2.3", "res": [-2.3, "1", "001", "31c"]},
        {"str": "1001.31c\t-2.3", "res": [-2.3, "1", "001", "31c"]},
        {"str": "15205 1", "res": [1, "15", "205", None]},
    ]

    def test_fromstring(self):
        """
        Test the creation of zaids from strings

        Returns
        -------
        None.

        """

        for test in self.tests:
            text = test["str"]
            zaid = Zaid.from_string(text)
            res = test["res"]
            assert zaid.fraction == res[0]
            assert zaid.element == res[1]
            assert zaid.isotope == res[2]
            assert zaid.library == res[3]

    def test_to_text(self):
        zaid = Zaid(
            "-2.3",
            "1",
            "001",
            "31c",
        )
        assert (
            zaid.to_text().strip()
            == "1001.31c       -2.300000E+0     $        WEIGHT(%)  AB(%)"
        )
        zaid.elem_mass_fraction = 0.1116316166
        zaid.ab = 0.514343484384

        assert (
            zaid.to_text().strip()
            == "1001.31c       -2.300000E+0     $        WEIGHT(%) 11.163 AB(%) 0.51434"
        )


class TestElement:
    zaid_strings = ["1001.31c   -1", "1002.31c   -3"]

    def _buildElem(self):
        zaids = []
        for zaidstr in self.zaid_strings:
            zaids.append(Zaid.from_string(zaidstr))

        elem = Element(zaids)
        return elem

    def test_update_zaidinfo(self):
        """
        Test ability to get additional info for the zaids
        """

        elem = self._buildElem()

        # Check for the correct element
        elem.Z = "1"

        # Check the correct update of infos in element
        elem.update_zaidinfo(LIBMAN, 0.5)
        res = [{"fullname": "H-1", "ab": 25}, {"fullname": "H-2", "ab": 75}]
        for zaid, checks in zip(elem.zaids, res):
            assert int(zaid.ab) == checks["ab"]
            assert zaid.fullname == checks["fullname"]

    def test_get_fraction(self):
        """
        Test correct recovery of element fraction
        """

        elem = self._buildElem()
        assert elem.get_fraction() == -4


class TestSubmaterial:
    def test_get_info(self):
        txt = [
            "C header",
            "8016.31c        1.333870E-2     $ O-16   AB(%) 99.757",
            "8017.31c        5.081060E-6     $ O-17   AB(%) 0.038",
            "8018.31c        2.741100E-5     $ O-18   AB(%) 0.205",
        ]

        submat = SubMaterial.from_text(txt)
        df_el, df_zaid = submat.get_info(LIBMAN)
        assert len(df_el) == 1
        assert len(df_zaid) == 3
        assert len(df_el.columns) == 2
        assert len(df_zaid.columns) == 3

    def test_dosimetry(self):
        lm = LibManager(dosimetry_lib="21c")
        txt = [
            "C header",
            "8016.31c        1     $ O-16   AB(%) 99.757",
            "8017.31c        1     $ O-17   AB(%) 0.038",
            "8018.21c        1     $ O-18   AB(%) 0.205",
        ]
        material = SubMaterial.from_text(txt)
        material.translate("00c", lm)
        assert material.zaidList[2].library == "21c"
        assert material.zaidList[0].library == "00c"

    def test_zaid_not_found(self):
        txt = [
            "C header",
            "8001.31c        1     $ O-16   AB(%) 99.757",
        ]
        material = SubMaterial.from_text(txt)
        material.translate("00c", LIBMAN)


class TestMaterial:
    def test_natural_expansion(self):
        lm = LibManager()
        material = Material.from_text(["C Header", "M1 1000.31c   -2.3"])
        material.translate("31c", lm)
        assert len(material.submaterials[0].zaidList) == 2

    def test_switch_fraction(self):
        """test needs to be conducted both for the creation through
        from_input() and from the creation from the input class
        """
        # Read a material
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)

        for matcard in [matcard1, matcard2]:
            # Fake translation in order to normalize the fractions
            material = matcard[0]
            # material._update_info(LIBMAN)
            original = material.to_text()

            # -- Switch back and forth --
            # this first one should do nothing
            material.switch_fraction("atom", LIBMAN)
            unchanged = material.to_text()
            assert original == unchanged
            # switch to mass
            material.switch_fraction("mass", LIBMAN)
            mass = material.to_text()
            # verify that they are indeed different
            assert original != mass
            # change again, should do nothing
            material.switch_fraction("mass", LIBMAN)
            unchanged = material.to_text()
            assert unchanged == mass
            # go back to atom
            material.switch_fraction("atom", LIBMAN)
            atom = material.to_text()
            # at last check that the inplace option works
            material.switch_fraction("mass", LIBMAN, inplace=False)
            inplace = material.to_text()
            assert inplace == atom
            # go back to mass
            material.switch_fraction("mass", LIBMAN)
            massnorm = material.to_text()
            assert massnorm == mass

    def test_switch_pnnl(self):
        # --- Test the PNNL with Bismuth Germanate (BGO) ---
        # read the material cards
        with as_file(resources.joinpath("BGO_mass.i")) as inp:
            matcard = MatCardsList.from_input(inp)
        mass_material = matcard[0]

        with as_file(resources.joinpath("BGO_atom.i")) as inp:
            matcard = MatCardsList.from_input(inp)
        atom_material = matcard[0]

        # Switch the mass fraction to atomic fraction
        mass_material.switch_fraction("atom", LIBMAN)
        print(mass_material.to_text())

        tolerance = 1e-5  # tolerance for the difference with respect to pnnl
        switched_sub = mass_material.submaterials[0]
        pnnl_sub = atom_material.submaterials[0]

        for zaid1, zaid2 in zip(switched_sub.zaidList, pnnl_sub.zaidList):
            diff = zaid1.fraction - zaid2.fraction
            assert diff < tolerance

    def test_from_zaids(self):
        zaids = [("1001", -100), ("B-0", -200), ("C-12", -50)]
        mat = Material.from_zaids(zaids, LIBMAN, "31c", "header")
        assert len(mat.submaterials[0].zaidList) == 4
        zaid = mat.submaterials[0].zaidList[0]
        assert zaid.element == "1"
        assert zaid.isotope == "001"
        assert zaid.fraction == -100.0

        zaid = mat.submaterials[0].zaidList[1]
        assert zaid.element == "5"
        assert zaid.isotope == "010"

        zaid = mat.submaterials[0].zaidList[2]
        assert zaid.element == "5"
        assert zaid.isotope == "011"

        zaids = [
            (1000, -4.7),
            (5000, -30.4),
            (6000, -28.3),
            (11000, -3.2),
            (16000, -33.1),
            (14000, -0.06),
            (26000, -0.08),
            (7000, -0.4),
        ]
        mat = Material.from_zaids(zaids, LIBMAN, "31c", "header")
        txt = """C header
M1
       1001.31c       -4.698638E+0     $ H-1    AB(%) 99.971
       1002.31c       -1.361756E-3     $ H-2    AB(%) 0.028974
       5010.31c       -5.531343E+0     $ B-10   AB(%) 18.195
       5011.31c       -2.486866E+1     $ B-11   AB(%) 81.805
       6012.31c       -2.797523E+1     $ C-12   AB(%) 98.852
       6013.31c       -3.247744E-1     $ C-13   AB(%) 1.1476
      11023.31c       -3.200000E+0     $ Na-23  AB(%) 100.0
      16032.31c       -3.130391E+1     $ S-32   AB(%) 94.574
      16033.31c       -2.596888E-1     $ S-33   AB(%) 0.78456
      16034.31c       -1.530534E+0     $ S-34   AB(%) 4.624
      16036.31c       -5.866145E-3     $ S-36   AB(%) 0.017722
      14028.31c       -5.513970E-2     $ Si-28  AB(%) 91.899
      14029.31c       -2.892181E-3     $ Si-29  AB(%) 4.8203
      14030.31c       -1.968119E-3     $ Si-30  AB(%) 3.2802
      26054.31c       -4.516447E-3     $ Fe-54  AB(%) 5.6456
      26056.31c       -7.352122E-2     $ Fe-56  AB(%) 91.902
      26057.31c       -1.728295E-3     $ Fe-57  AB(%) 2.1604
      26058.31c       -2.340355E-4     $ Fe-58  AB(%) 0.29254
       7014.31c       -3.983744E-1     $ N-14   AB(%) 99.594
       7015.31c       -1.625644E-3     $ N-15   AB(%) 0.40641"""
        compare_without_dollar_comments(txt, mat.to_text())

    def test_get_tad(self):
        with as_file(resources.joinpath("tad_test.i")) as inp:
            mcnp_inp = Input.from_input(inp)
            materials = mcnp_inp.materials
        dens_mats = {
            "316L(N)-IG": 7.93,
            "CuCrZr": 8.9,
        }
        mats = {
            "316L(N)-IG": 100,
            "CuCrZr": 302,
        }
        tads = {
            "316L(N)-IG": 8.59386e-002,
            "CuCrZr": 8.46480e-002,
        }
        libman = LibManager()
        for mat, mat_num in mats.items():
            mat_card = materials["M" + str(mat_num)]
            assert (
                pytest.approx(mat_card.get_tad(dens_mats[mat], libman), 1e-4)
                == tads[mat]
            )

    def test_get_density(self):
        with as_file(resources.joinpath("tad_test.i")) as inp:
            mcnp_inp = Input.from_input(inp)
            materials = mcnp_inp.materials
        dens_mats = {
            "316L(N)-IG": 7.93,
            "CuCrZr": 8.9,
        }
        mats = {
            "316L(N)-IG": 100,
            "CuCrZr": 302,
        }
        tads = {
            "316L(N)-IG": 8.59386e-002,
            "CuCrZr": 8.46480e-002,
        }
        libman = LibManager()
        for mat, mat_num in mats.items():
            mat_card = materials["M" + str(mat_num)]
            assert (
                pytest.approx(mat_card.get_density(tads[mat], libman), 1e-4)
                == dens_mats[mat]
            )

    def test_get_info(self):
        txt = [
            "C header",
            "8016.31c        1.333870E-2     $ O-16   AB(%) 99.757",
            "8017.31c        5.081060E-6     $ O-17   AB(%) 0.038",
            "8018.31c        2.741100E-5     $ O-18   AB(%) 0.205",
        ]
        submat1 = SubMaterial.from_text(txt)
        txt = ["C header", "1001.31c  1e-2", "8016.31c  1e-2"]
        submat2 = SubMaterial.from_text(txt)
        material = Material(None, None, "M1", submaterials=[submat1, submat2])
        lm = LibManager()
        df, df_elem = material.get_info(lib_manager=lm, zaids=True)
        assert df.iloc[0]["Atom Fraction"] != df.iloc[0]["Mass Fraction"]

    def test_update_info(self):
        libman = LibManager()
        material = Material.from_zaids([(1000, 1), (56000, 1)], libman, "31c")
        assert (
            "1001.31c        9.998550E-1     $ H-1    WEIGHT(%) 0.72865 AB(%) 99.986"
            in material.to_text()
        )


class TestMatCardList:
    """test needs to be conducted both for the creation through
    from_input() and from the creation from the input class
    """

    def test_frominput(self):
        """
        Test basic properties

        Returns
        -------
        None.

        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)

        for matcard in [matcard1, matcard2]:
            assert len(matcard.materials) == 3
            assert len(matcard.matdic) == 3

    def test_headers(self):
        """
        test correct material headers reading

        Returns
        -------
        None.

        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)

        headers = {"m1": "C Header M1\n", "m2": "C Header M2\n", "m102": ""}
        for matcard in [matcard1, matcard2]:
            for key, header in headers.items():
                assert matcard[key].header == header

    def test_subheaders(self):
        """
        Test correct reading of submaterial headers

        Returns
        -------
        None.

        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)

        headers = {
            "m1": ["C M1-submat1", "C M1-Submat 2"],
            "m2": ["", "C M2-submat1\nC second line"],
            "m102": [""],
        }

        for matcard in [matcard1, matcard2]:
            for key, subheaders in headers.items():
                for i, submat in enumerate(matcard[key].submaterials):
                    assert submat.header == subheaders[i]

    def test_zaidnumbers(self):
        """
        Test correct number of zaids allocated in submaterials

        Returns
        -------
        None.

        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)

        zaids_dic = {"m1": [2, 1], "m2": [1, 1], "m102": [5]}

        for matcard in [matcard1, matcard2]:
            for key, zaids in zaids_dic.items():
                for i, submat in enumerate(matcard[key].submaterials):
                    assert len(submat.zaidList) == zaids[i]

    def test_translation(self):
        """
        Test that translation works (all possile modes)
        """
        # Dic mode 1
        matcard1 = deepcopy(activation_matcard1)
        matcard2 = deepcopy(activation_matcard2)

        newlib = {"21c": "31c", "99c": "81c"}
        for matcard in [matcard1, matcard2]:
            matcard.translate(newlib, LIBMAN)
            translation = matcard.to_text()
            assert translation.count("31c") == 3
            assert translation.count("81c") == 3

        # dic mode 2 - test 1
        matcard1 = deepcopy(activation_matcard1)
        matcard2 = deepcopy(activation_matcard2)
        newlib = {"99c": ["1001"], "21c": ["28061", "28062", "28064", "29063", "5010"]}
        for matcard in [matcard1, matcard2]:
            matcard.translate(newlib, LIBMAN)
            translation = matcard.to_text()
            assert translation.count("99c") == 0
            assert translation.count("21c") == 5
            assert translation.count("81c") == 1

        # dic mode 2 - test 2
        matcard1 = deepcopy(activation_matcard1)
        matcard2 = deepcopy(activation_matcard2)
        newlib = {"99c": ["1001"], "21c": ["28061", "28062", "28064", "29063"]}
        for matcard in [matcard1, matcard2]:
            try:
                matcard.translate(newlib, LIBMAN)
                assert False
            except ValueError:
                assert True

        # classic mode
        matcard1 = deepcopy(inp2_matcard1)
        matcard2 = deepcopy(inp2_matcard2)
        for matcard in [matcard1, matcard2]:
            matcard.translate("21c", LIBMAN)
            translation = matcard.to_text()
            assert translation.count("21c") == 10

    def test_get_info(self):
        """
        Barely tests that everything is created
        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)
        for matcard in [matcard1, matcard2]:
            df, df_elem = matcard.get_info(LIBMAN, zaids=True)
            assert len(df) == 8
            assert len(df_elem) == 7

    def test_generate_material(self):
        # using atom fraction
        matcard = deepcopy(inp3_matcard2)
        materials = ["m1", "M2"]
        percentages = [0.5, 0.5]
        newlib = "31c"
        # using atom fraction
        fraction_type = "atom"
        newmat = matcard.generate_material(
            materials, percentages, newlib, LIBMAN, fractiontype=fraction_type
        )
        fileA = os.path.join(resources, "newmat_atom")
        text_A = ""
        with open(fileA, "r") as infile:
            for line in infile:
                text_A = text_A + line

        compare_without_dollar_comments(text_A, newmat.to_text())

        # using mass fraction
        fraction_type = "mass"
        newmat = matcard.generate_material(
            materials, percentages, newlib, LIBMAN, fractiontype=fraction_type
        )
        fileB = os.path.join(resources, "newmat_mass")
        text_B = ""
        with open(fileB, "r") as infile:
            for line in infile:
                text_B = text_B + line

        compare_without_dollar_comments(text_B, newmat.to_text())

    def test_fractions_to_atom_density(self):
        matcard = deepcopy(inp3_matcard2)
        density = 2

        submats = matcard["m1"].switch_fraction("mass", LIBMAN, inplace=False)
        newmat = Material(None, None, "", submaterials=submats)
        submats = newmat.switch_fraction("atom", LIBMAN, inplace=False)
        fraction = submats[0].zaidList[0].fraction

        tad = matcard["m1"].get_tad(density, LIBMAN)

        matcard.fractions_to_atom_densities(LIBMAN, density)

        assert (
            pytest.approx(matcard["m1"].submaterials[0].zaidList[0].fraction)
            == fraction * tad
        )
        assert matcard["m1"].submaterials[0].zaidList[1].fraction > 0

    def test_generate_material_from_created_materials(self):
        materials = MatCardsList([])
        zaids = [("1001", -100), ("B-0", -200), ("C-12", -50)]
        newmat = Material.from_zaids(zaids, LIBMAN, "31c", "B4CH2", 302)
        materials.append(newmat)
        zaids = [("1001", 2), ("O-0", 1)]
        newmat = Material.from_zaids(zaids, LIBMAN, "31c", "Water", 400)
        materials.append(newmat)
        newmat = materials.generate_material(
            ["M302", "M400"],
            [0.2, 0.8],
            "31c",
            LIBMAN,
            fractiontype="atom",
            mat_name="M500",
        )
        assert (
            newmat.to_text()
            == """C Material: M302 Percentage: 20.0% (atom)
C Material: M400 Percentage: 80.0% (atom)
M500
C B4CH2
C M302, submaterial 1
C no submat header
       1001.31c        1.628143E-1     $ H-1    WEIGHT(%) 28.571 AB(%) 100.0
       5010.31c        5.963525E-3     $ B-10   WEIGHT(%) 57.143 AB(%) 19.65
       5011.31c        2.438520E-2     $ B-11   WEIGHT(%) 57.143 AB(%) 80.35
       6012.31c        6.837012E-3     $ C-12   WEIGHT(%) 14.286 AB(%) 100.0
C Water
C M400, submaterial 1
C no submat header
       1001.31c        5.333336E-1     $ H-1    WEIGHT(%) 11.189 AB(%) 100.0
       8016.31c        2.660188E-1     $ O-16   WEIGHT(%) 88.811 AB(%) 99.757
       8017.31c        1.022667E-4     $ O-17   WEIGHT(%) 88.811 AB(%) 0.03835
       8018.31c        5.453336E-4     $ O-18   WEIGHT(%) 88.811 AB(%) 0.2045"""
        )


def compare_without_dollar_comments(text_A: str, text_B: str):
    # strip all dollars
    pattern = re.compile(r"\$.*?\n")
    newtextA = pattern.sub("", text_A)
    newtextA = newtextA.split("$")[0]
    newtextB = pattern.sub("", text_B)
    newtextB = newtextB.split("$")[0]
    assert newtextA == newtextB
