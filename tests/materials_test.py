import os
from copy import deepcopy
from importlib.resources import files, as_file

from f4eparser.input.materials import (Element, Zaid, MatCardsList)
from f4eparser.input.libmanager import LibManager
from f4eparser.input.inputAPI import Input
import f4eparser.resources as pkg_res
import tests.resources.materials as mat_res


resources = files(mat_res)
resources_pkg = files(pkg_res)

XSDIR = as_file(resources.joinpath('xsdir_mcnp6.2'))
ISOTOPES_FILE = as_file(resources_pkg.joinpath('Isotopes.txt'))
# Other
with (XSDIR as xsdir_file,
      ISOTOPES_FILE as isotopes_file):

    LIBMAN = LibManager(xsdir_file,
                        isotopes_file=isotopes_file, defaultlib='81c')

# Files
with as_file(resources.joinpath('mat_test.i')) as inp:
    inp_matcard1 = MatCardsList.from_input(inp)
    inp_matcard2 = Input.from_input(inp).materials

with as_file(resources.joinpath('mat_test2.i')) as inp2:
    inp2_matcard1 = MatCardsList.from_input(inp2)
    inp2_matcard2 = Input.from_input(inp2).materials

with as_file(resources.joinpath('test.i')) as inp3:
    inp3_matcard1 = MatCardsList.from_input(inp3)
    inp3_matcard2 = Input.from_input(inp3).materials

with as_file(resources.joinpath('activation.i')) as inp:
    activation_matcard1 = MatCardsList.from_input(inp)
    activation_matcard2 = Input.from_input(inp).materials


class TestZaid:

    tests = [{'str': '1001.31c   -2.3', 'res': [-2.3, '1', '001', '31c']},
             {'str': '1001.31c\t-2.3', 'res': [-2.3, '1', '001', '31c']},
             {'str': '15205 1', 'res': [1, '15', '205', None]}]

    def test_fromstring(self):
        """
        Test the creation of zaids from strings

        Returns
        -------
        None.

        """

        for test in self.tests:
            text = test['str']
            zaid = Zaid.from_string(text)
            res = test['res']
            assert zaid.fraction == res[0]
            assert zaid.element == res[1]
            assert zaid.isotope == res[2]
            assert zaid.library == res[3]


class TestElement:
    zaid_strings = ['1001.31c   -1', '1002.31c   -3']

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
        elem.Z = '1'

        # Check the correct update of infos in element
        elem.update_zaidinfo(LIBMAN)
        res = [{'fullname': 'H-1', 'ab': 25},
               {'fullname': 'H-2', 'ab': 75}]
        for zaid, checks in zip(elem.zaids, res):
            assert int(zaid.ab) == checks['ab']
            assert zaid.fullname == checks['fullname']

    def test_get_fraction(self):
        """
        Test correct recovery of element fraction
        """

        elem = self._buildElem()
        assert elem.get_fraction() == -4


class TestMaterial:

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
            material.update_info(LIBMAN)
            original = material.to_text()

            # -- Switch back and forth --
            # this first one should do nothing
            material.switch_fraction('atom', LIBMAN)
            unchanged = material.to_text()
            assert original == unchanged
            # switch to mass
            material.switch_fraction('mass', LIBMAN)
            mass = material.to_text()
            # change again, should do nothing
            material.switch_fraction('mass', LIBMAN)
            unchanged = material.to_text()
            assert unchanged == mass
            # go back to atom
            material.switch_fraction('atom', LIBMAN)
            atom = material.to_text()
            # at last check that the inplace oprion works
            material.switch_fraction('mass', LIBMAN, inplace=False)
            inplace = material.to_text()
            assert inplace == atom
            # go back to mass
            material.switch_fraction('mass', LIBMAN)
            massnorm = material.to_text()
            assert massnorm == mass

    def test_switch_pnnl(self):
        # --- Test the PNNL with Bismuth Germanate (BGO) ---
        # read the material cards
        with as_file(resources.joinpath('BGO_mass.i')) as inp:
            matcard = MatCardsList.from_input(inp)
        mass_material = matcard[0]

        with as_file(resources.joinpath('BGO_atom.i')) as inp:
            matcard = MatCardsList.from_input(inp)
        atom_material = matcard[0]

        # Switch the mass fraction to atomic fraction
        mass_material.switch_fraction('atom', LIBMAN)
        print(mass_material.to_text())

        tolerance = 1e-5  # tolerance for the difference with respect to pnnl
        switched_sub = mass_material.submaterials[0]
        pnnl_sub = atom_material.submaterials[0]

        for zaid1, zaid2 in zip(switched_sub.zaidList, pnnl_sub.zaidList):
            diff = zaid1.fraction - zaid2.fraction
            assert diff < tolerance


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

        headers = {'m1': 'C Header M1\n', 'm2': 'C Header M2\n', 'm102': ''}
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

        headers = {'m1': ['C M1-submat1', 'C M1-Submat 2'],
                   'm2': ['', 'C M2-submat1\nC second line'],
                   'm102': ['']}

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

        zaids_dic = {'m1': [2, 1],
                     'm2': [1, 1],
                     'm102': [5]}

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

        newlib = {'21c': '31c', '99c': '81c'}
        for matcard in [matcard1, matcard2]:
            matcard.translate(newlib, LIBMAN)
            translation = matcard.to_text()
            assert translation.count('31c') == 3
            assert translation.count('81c') == 3

        # dic mode 2 - test 1
        matcard1 = deepcopy(activation_matcard1)
        matcard2 = deepcopy(activation_matcard2)
        newlib = {'99c': ['1001'], '21c': ['28061', '28062', '28064', '29063',
                                           '5010']}
        for matcard in [matcard1, matcard2]:
            matcard.translate(newlib, LIBMAN)
            translation = matcard.to_text()
            assert translation.count('99c') == 0
            assert translation.count('21c') == 5
            assert translation.count('81c') == 1

        # dic mode 2 - test 2
        matcard1 = deepcopy(activation_matcard1)
        matcard2 = deepcopy(activation_matcard2)
        newlib = {'99c': ['1001'], '21c': ['28061', '28062', '28064', '29063']}
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
            matcard.translate('21c', LIBMAN)
            translation = matcard.to_text()
            assert translation.count('21c') == 10

    def test_get_info(self):
        """
        Barely tests that everything is created
        """
        matcard1 = deepcopy(inp_matcard1)
        matcard2 = deepcopy(inp_matcard2)
        for matcard in [matcard1, matcard2]:
            df, df_elem = matcard.get_info(LIBMAN, zaids=True, complete=True)
            assert len(df) == 8
            assert len(df_elem) == 7

    def test_generate_material(self):
        # using atom fraction
        matcard = deepcopy(inp3_matcard2)
        materials = ['m1', 'M2']
        percentages = [0.5, 0.5]
        newlib = '31c'
        # using atom fraction
        fraction_type = 'atom'
        newmat = matcard.generate_material(materials, percentages, newlib,
                                           LIBMAN,
                                           fractiontype=fraction_type)
        fileA = os.path.join(resources, 'newmat_atom')
        text_A = ''
        with open(fileA, 'r') as infile:
            for line in infile:
                text_A = text_A+line

        assert text_A == newmat.to_text()

        # using mass fraction
        fraction_type = 'mass'
        newmat = matcard.generate_material(materials, percentages, newlib,
                                           LIBMAN,
                                           fractiontype=fraction_type)
        fileB = os.path.join(resources, 'newmat_mass')
        text_B = ''
        with open(fileB, 'r') as infile:
            for line in infile:
                text_B = text_B+line

        assert text_B == newmat.to_text()
