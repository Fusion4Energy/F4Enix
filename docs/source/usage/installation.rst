.. _install:

############
Installation
############

.. warning::
    Please remeber that python versions lower than **v3.10** are not supported.

Both Windows and Linux OS are supported.


Install from PyPI
=================
*Coming soon*


Local install
=============
To have a simple local installation (default choice).

Navigate to the F4Eparser folder and type the following to be sure that you
pip is up to date:

``python -m pip install --upgrade pip setuptools wheel``

then run in your python environment of choice:

``python -m pip install .`` [Windows]

``python3 -m pip install .`` [Linux]

Developer mode install
======================

This installation will only save the package info in your local environment
while it will point to your local folder to get the actual package code.
This means that when you have the package installed (i.e. that you can
use in whatever directory) you will still be able to modify the package
source code without the need for re-installing the package.

A common use for this is to clone the F4Enix GitHub repository on your local
pc and then install it from there. In this way, it will be possible to edit
the code (either by manual modification or pulling from online repo) without
having to reinstall the package. 

This also means though that if you move the F4Enix folder, the link will
be broken and you will need to reinstall the package.

Navigate to the F4Eparser folder and type the following to be sure that you
pip is up to date:

``python -m pip install --upgrade pip setuptools wheel``

and then to install the package:

``python -m pip install -e .`` [Windows]

``python3 -m pip install -e .`` [Linux]


.. important:: 
    if you are encountering any sort of trouble, the first action when
    troubleshooting would be to create a new fresh virtual environment in order
    to be sure to avoid clashes between previously installed python packages.
    In conda this would be done with:
    ``conda create -n <env_name> python=3.10``