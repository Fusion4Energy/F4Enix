[![Testing windows](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_win.yml/badge.svg?branch=main)](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_win.yml)
[![Testing linux](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_linux.yml/badge.svg?branch=main)](https://github.com/Radiation-Transport/F4Enix/actions/workflows/AutomatedTests_linux.yml)

# F4Enix
Parser for Monte Carlo simulations input and output files

**Pyhton >3.10!**

Both Windows and Linux supported.

## Install
### Disclaimer
Currently the package can be installed only in developer mode.

The step 0 of this procedure, that is, if you are encountering any
sort of trouble, would be to create a new fresh virtual environment.
In conda this would be done with:
```
conda create -n <env_name> python=3.10
```
Please remeber that python versions lower than 3.10 are not supported.

### Proper installation
Since the wheels of the package are not being built, the package can be
install in editing (developer) mode. In order to do so, navigate to the
F4Eparser folder and type the following to be sure that you pip is
up to date:
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