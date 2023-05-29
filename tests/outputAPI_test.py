from importlib.resources import files, as_file

from f4enix.output.outputAPI import Output
import tests.resources.output as resources

RESOURCES = files(resources)


class TestOutput:
    def test_print_lp_debug(self, tmpdir):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)

        # This has no LP
        outp.print_lp_debug(tmpdir)
        checks = outp.get_statistical_checks()
        assert len(checks) == 12

        with as_file(RESOURCES.joinpath('out_lp.txt')) as file:
            outp = Output(file)
        # This has LP
        outp.print_lp_debug(tmpdir)

    def test_get_NPS(self):
        with as_file(RESOURCES.joinpath('test_o')) as file:
            outp = Output(file)
        assert outp.get_NPS() == int(1e4)
