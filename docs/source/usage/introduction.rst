What is F4Enix
==============

The general idea, is to develop a series of specific python modules
responsible for the parsing and manipulation of all the main types of files
that are commonly produced during a complex Monte Carlo particles and
radiation transport analysis. These are divided in two main subpackages which
are ``f4enix.input`` and ``f4enix.output``.

The following is a list of the kind of files that are currently at least
partially supported by F4Enix grouped by codes.

**MCNP:**

* MCNP Input file, which essentially is a wrapper of some parts of the `numjuggler <https://numjuggler.readthedocs.io/>`_ python package.
  (`Examples <examples/input/jupyters/mcnp_inp.html#mcnp-input-files>`_)
* MCNP Output file (`Examples <examples/output/jupyters/mcnp_out.html#mcnp-output-files>`_)
* MCTAL file (`Examples <examples/output/jupyters/mctal.html#mctal-file>`_)
* MESHTAL file, including modified versions produced by D1SUNED
  (`Examples <examples/output/jupyters/meshtal.html#meshtal-file>`_)
* RSSA file (`Examples <examples/output/jupyters/rssa.html#rssa-file>`_)
* EEOUT file (`Examples <examples/output/jupyters/eeout.html#eeout-files-unstructured-meshes>`_)
* WWINP Weight-Windows file (`Examples <examples/input/jupyters/ww.html#example-of-use-via-python-scripting>`_)
* GVR capabilities (`Examples <examples/input/jupyters/ww.html#creation-of-a-gvr>`_)

**D1SUNED:**

* D1SUNED input file (`Examples <examples/input/jupyters/d1suned.html#d1s-uned-files>`_)
* D1SUNED meshinfo file (`Examples <examples/output/jupyters/meshinfo.html#meshinfo-files>`_)

**FISPACT:**

* FISPACT legacy output (for pathways) (`Examples <examples/output/jupyters/fispact_legacy.html#fispact-ii-legacy-output>`_)
* Important Pathways library for ITER applications (`Examples <examples/output/jupyters/pathwaylib.html#decay-pathway-library>`_)


All classes and methods of the F4Enix API are documented and usage examples
are provided in all the most important classes documentations. Additionally more
structured examples of pre and post-processing pipelines are provided in the form
of compiled jupyter notebooks.
Everything that involves mesh output is dealt with the very versatile python
package `PyVista <https://docs.pyvista.org/version/stable/index.html>`_.

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

Python has been chosen as the programming language for making it easy to use,
portable, and easy to be intagrated with other scientific libraries and tools.

Continous Integration (CI) procedures are implemented for F4Enix, additional
details may be found at :ref:`CI details`.
Both Linux and Windows OS are supported and tested.

**Short/mid term goals for the project:**

* general issues fixing
* improve documentation

**Long term goals for the project:**

* compatibilty with newer output formats brought by MCNP v6.3
* increase features in all modules depending on needs
* change MCNP input parser engine from ``numjuggler`` to something more robust
