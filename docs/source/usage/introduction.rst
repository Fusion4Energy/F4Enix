############
Introduction
############

Why F4Enix
==========

The reliable computation of nuclear responses for ITER-like systems and
components is a complex and resource-intensive process.
It encompasses the preparation of exceptionally large and detailed computer
models and the use of sophisticated, non-commercial scientific software
(e.g. MCNP) in uncommon large-scale parallel computing facilities (HPC).

To address the challenges associated with this process, the F4Enix python
package has been developed by the neutronics team of Fusion For Energy (F4E).
The primary objective of F4Enix is to automate and streamline the pre and 
post-processing operations involved in nuclear response computations for
ITER or similar projects. By leveraging the power of Python, this package
provides a comprehensive set of tools for efficient and high-quality parsing
and manipulation of MCNP inputs and outputs. These tools aim to significantly
enhance the efficiency, capability, and overall quality of the entire nuclear
analysis workflow.

To foster collaboration, encourage improvement, and avoid duplication of
efforts, the development of F4Enix follows an open-source approach.
The open-source nature of the project ensures accessibility at a pan-European
level and facilitates engagement with a wide user community. It also enables
users to contribute to the debugging and enhancement of the package,
ensuring continuous development and improvement.

The package is made available on GitHub at [FUTURE LINK].

Python has been chosen as the programming language for making it easy to use,
portable, and easy to be intagrated with other scientific libraries and tools.
 

What is F4Enix
==============

The general idea, is to develop a series of specific python modules
responsible for the parsing and manipulation of all the main types of files
that are commonly produced during a complex Monte Carlo particles and 
radiation transport analysis. These are divided in two main subpackages which
are ``f4enix.input`` and ``f4enix.output``.

The following is a list of the kind of files that are currently at least partially
supported by F4Enix:

* MCNP input file by :py:mod:`f4enix.input.inputAPI`
* MCNP MCTAL file by :py:mod:`f4enix.output.mctal`
* MCNP MESHTAL file, including modified versions produced by D1SUNED
  by :py:mod:`f4enix.output.meshtal`
* MCNP output file (only some basic features) by :py:mod:`f4enix.output.outputAPI`
* VTK files by :py:mod:`f4enix.output.pyvistawrap`

Short/mid term goals for the project:

* support for MCNP Weight-Windows files
* support for MCNP unstructured meshes output
* support for RSSA files
* support for the meshinfo file from CuV FMESH approach
* build a common plotter to be potentially used by different modules
* increase automatic test coverage
* improve documentation