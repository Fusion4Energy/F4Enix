import sys
import os
import pytest
from importlib.resources import files, as_file

from f4enix.input.d1suned import (Reaction, ReactionFile, Irradiation,
                                  IrradiationFile, REACFORMAT)
from f4enix.input.libmanager import LibManager
import tests.resources.d1suned as res

RESOURCES = files(res)
# INP = os.path.join(cp, 'TestFiles', 'parserD1S', 'reac_fe')
# Files


class TestIrradiationFile:

    @pytest.mark.parametrize('file', ['irr_test', 'irr_test2'])
    def test_fromtext(self, file):
        """
        Test parsing irradiation file
        """
        with as_file(RESOURCES.joinpath(file)) as inp:
            irrfile = IrradiationFile.from_text(inp)
        if file == 'irr_test':
            self._assert_file1(irrfile)
        elif file == 'irr_test2':
            self._assert_file2(irrfile)

    @staticmethod
    def _assert_file1(irrfile: IrradiationFile):
        assert len(irrfile.irr_schedules) == 4
        TestIrradiation.assert_irr(irrfile.irr_schedules[0])

    @staticmethod
    def _assert_file2(irrfile: IrradiationFile):
        assert len(irrfile.irr_schedules) == 4
        TestIrradiation.assert_irr(irrfile.irr_schedules[0])

    @pytest.mark.parametrize('file', ['irr_test', 'irr_test2'])
    def test_write(self, tmpdir, file):
        """
        Test writing irradiation file 1
        """
        with as_file(RESOURCES.joinpath(file)) as inp:
            irrfile = IrradiationFile.from_text(inp)
        irrfile.write(tmpdir)
        outfile = os.path.join(tmpdir, irrfile.name)
        irrfile = IrradiationFile.from_text(outfile)
        if file == 'irr_test':
            self._assert_file1(irrfile)
        elif file == 'irr_test2':
            self._assert_file2(irrfile)

    def test_get_daughters(self):
        with as_file(RESOURCES.joinpath('irr_test')) as inp:
            irrfile = IrradiationFile.from_text(inp)
        daughters = irrfile.get_daughters()
        assert daughters == ['24051', '25054', '26055', '26059']

    def test_get_irrad(self):
        with as_file(RESOURCES.joinpath('irr_test')) as inp:
            irrfile = IrradiationFile.from_text(inp)
        # Check the None
        irradiation = irrfile.get_irrad('20051')
        assert irradiation is None
        # Check true value
        irradiation = irrfile.get_irrad('26055')
        assert irradiation.daughter == '26055'


class TestIrradiation:

    def test_reading(self):
        """
        Test the reading of irradiation line
        """
        text = '   24051     2.896e-07    5.982e+00    5.697e+00     Cr51'
        irr = Irradiation.from_text(text, 2)
        self.assert_irr(irr)

    @staticmethod
    def assert_irr(irr: Irradiation):
        """
        Assert irradiation
        """
        assert irr.daughter == '24051'
        assert irr.lambd == '2.896e-07'
        assert irr.times[0] == '5.982e+00'
        assert irr.times[1] == '5.697e+00'
        assert irr.comment == 'Cr51'

    def test_equivalence(self):
        # Equivalent
        text = '   24051     2.896e-07    5.982e+00    5.697e+00     Cr51'
        irr1 = Irradiation.from_text(text, 2)
        text = '   24051     2.896e-07    5.982e+00    5.697     '
        irr2 = Irradiation.from_text(text, 2)
        assert irr1 == irr2

        # Not equal
        text = '   24051     2.896e-07    5.697e+00    5.982e+00     Cr51'
        irr3 = Irradiation.from_text(text, 2)
        text = '   24051     2.896e-07    5.697e+00    Cr51'
        irr4 = Irradiation.from_text(text, 1)
        assert irr1 != irr3
        assert irr1 != {}
        assert irr1 != irr4


class TestReaction:

    def test_fromtext1(self):
        """
        Test different formatting possibilities
        """
        text = '   26054.99c  102  26055     Fe55'
        reaction = Reaction.from_text(text)
        assert reaction.parent == '26054.99c'
        assert reaction.MT == '102'
        assert reaction.daughter == '26055'
        assert reaction.comment == 'Fe55'

    def test_fromtext2(self):
        """
        Test different formatting possibilities
        """
        text = '26054.99c 102   26055 Fe55  and some'
        reaction = Reaction.from_text(text)
        assert reaction.parent == '26054.99c'
        assert reaction.MT == '102'
        assert reaction.daughter == '26055'
        assert reaction.comment == 'Fe55 and some'

    def test_changelib(self):
        """
        Test change library tag
        """
        rec = Reaction('26054.99c', '102', '26055')
        rec.change_lib('31c')
        assert rec.parent == '26054.31c'

    def test_write(self):
        """
        check writing
        """
        text = '26054.99c  102  26055     Fe55 and  some'
        reaction = Reaction.from_text(text)
        ftext = reaction._get_text()
        comptext = ['26054.99c', '102', '26055', 'Fe55 and some']
        assert comptext == ftext


class TestReactionFile:

    @pytest.fixture
    def lm(self):
        # xsdirpath = os.path.join(cp, 'TestFiles', 'libmanager', 'xsdir')
        # isotopes_file = os.path.join(root, 'jade', 'resources', 'Isotopes.txt')
        return LibManager()

    def test_fromtext(self):
        """
        right number of reactions
        """
        with as_file(RESOURCES.joinpath('reac_fe')) as inp:
            reac_file = ReactionFile.from_text(inp)
        print(reac_file.reactions)
        assert len(reac_file.reactions) == 10

    def test_write(self, tmpdir):
        """
        writing works
        """
        with as_file(RESOURCES.joinpath('reac_fe')) as inp:
            reac_file = ReactionFile.from_text(inp)

        reac_file.write(tmpdir)
        outpath = os.path.join(tmpdir, 'react')
        newfile = ReactionFile.from_text(outpath)
        # Remove the temporary file
        os.remove(outpath)
        # do some operations
        newfile.change_lib('31c')
        assert len(newfile.reactions) == 10
        # Check also first line
        rx = newfile.reactions[0]
        assert rx.parent == '26054.31c'
        assert rx.MT == '102'
        assert rx.daughter == '26055'
        assert rx.comment == 'Fe55'

    def test_translation(self, lm: LibManager):
        """
        test translation with libmanager where parents are available

        """
        newlib = '98c'

        with as_file(RESOURCES.joinpath('reac_fe')) as inp:
            reac_file = ReactionFile.from_text(inp)

        reac_file.change_lib(newlib, libmanager=lm)

        for reaction in reac_file.reactions:
            assert reaction.parent[-3:] == newlib

    def test_translation2(self, lm: LibManager):
        """
        test translation with libmanager where parents are not available

        """

        with as_file(RESOURCES.joinpath('reac2')) as inp:
            reac_file = ReactionFile.from_text(inp)

        newlib = '99c'

        reac_file.change_lib(newlib, libmanager=lm)

        for reaction in reac_file.reactions:
            assert reaction.parent[-3:] != newlib

    def test_get_parents(self):
        with as_file(RESOURCES.joinpath('reac_fe')) as inp:
            reac_file = ReactionFile.from_text(inp)
        parents = reac_file.get_parents()
        assert parents == ['26054', '26056', '26057', '26058']
