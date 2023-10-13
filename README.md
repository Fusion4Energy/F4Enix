[![Testing windows](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_win.yml/badge.svg?branch=main)](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_win.yml)
[![Testing linux](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_linux.yml/badge.svg?branch=main)](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_linux.yml)
[![PyPi version](https://badgen.net/pypi/v/f4enix/)](https://pypi.org/project/f4enix)
[![Documentation Status](https://readthedocs.org/projects/f4enix/badge/?version=latest)](https://f4enix.readthedocs.io/en/latest/?badge=latest)

# F4Enix
Parser for Monte Carlo simulations input and output files

**Pyhton >3.10!**

Both Windows and Linux supported.

Go to [F4Enix official documentation](https://f4enix.readthedocs.io/en/latest/) to get
more information on the library capabilities, examples and much more.

## Install
The easiest way to install F4Enix is using pip:

```
pip install f4enix
```

### Troubleshooting and developer mode
If any unexpected issue is encountered during installation, the first step for
its resolution would be to create a new fresh virtual environment.
In conda this would be done with:
```
conda create -n <env_name> python=3.10
```
Please remeber that python versions lower than 3.10 are not supported.

### Developer mode installation
There are many situations (especially while being an active developer) where it may be useful to install f4enix in developer mode. In order to do so, navigate to the
F4Eparser folder and type the following to be sure that your pip is
up to date (optional):
```
python -m pip install --upgrade pip setuptools wheel
```
and then to install the package:
```
python3 -m pip install -e .
```

This installation will only save the package info in your local environment
while it will point to your local folder to get the actual package code.
This means that when you have the package installed (i.e. that you can
use in whatever directory) you will still be able to modify the package
source code without the need for re-installing the package.

This also means though that if you move the F4Eparser folder, the link will
be broken and you will need to reinstall the package. Anyhow, if you are
working on your local GitHub repo to modify the code (as it should be) this
should never be moved anyway.