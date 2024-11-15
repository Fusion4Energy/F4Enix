from importlib.metadata import version

from .input.MCNPinput import Input
from .output.MCNPoutput import Output
from .output.mctal import Mctal
from .output.meshtal import Meshtal

__version__ = version("f4enix")
