######################
Contributing to F4Enix
######################

Thank you for considering contributing to F4Enix! This tool has been initially
created in order to better organize the different tools developed at F4E and
centralize its maintenance and documentation. F4E is a public entity and, as
such, whatever tool developed under its umbrella belongs to the public.
We look forward to community contributions in order for F4Enix to really be
a tool for the people by the people.

The purpose of this section is to document how the project is managed: how contributions
(bug fixes, enhancements, new features) are made, how they are evaluated,
who is permitted to merge pull requests, and what happens in the event of disagreements.
Once you have read through this section, the Development Workflow section outlines the
actual mechanics of making a contribution (forking, submitting a pull request, etc.).

These guidelines are inspired by what is done in the
`OpenMC <https://docs.openmc.org/en/stable/devguide/index.html>`_ community.

Terminology
===========
* A *Contributor* is any individual creating or commenting on an issue or pull request.
* A *Committer* is a subset of contributors who are authorized to review and merge pull requests.
* The *Project Lead Team* is a small group of committers who has the authority to mak final decisions
  and tag new releases.

Development Workflow
====================
Anyone wishing to make contributions to F4Enix should be fully acquainted and comfortable
working with git and GitHub. We assume here that you have git installed on your system,
have a GitHub account and that you are able to create/push to repositories on GitHub.
An easy way to approach the world of git actions and their integration with GitHub is to use
`GitHub Desktop <https://desktop.github.com/>`_.

Development of F4Enix relies heavily on branching; specifically, we use a branching model
sometimes referred to as git flow. If you plan to contribute to F4Enix development,
we highly recommend that you read this
`blog post <https://nvie.com/posts/a-successful-git-branching-model/>`_
to get a sense of how the branching
model works. There are two main branches that always exist: ``master`` and ``developing``.
The master branch is a stable branch that contains the latest release of the code.
The develop branch is where any ongoing development takes place prior to a release and is
not guaranteed to be stable. When the project lead team decides that a release should occur,
the ``developing`` branch is merged into master.

All new features, enhancements, and bug fixes should be developed on a branch that branches off
of develop. When the feature is completed, a pull request is initiated on GitHub that is
then reviewed by a committer. If the pull request is satisfactory, it is then merged into develop.
Note that a committer may not review their own pull request
(i.e., an independent code review is required).

Contribution steps
------------------
These steps apply to both new features and bug fixes. The general steps for contributing
are as follows:

#. Fork the main F4Enix repository from [TBD insert link here]. This will create a
   repository with the same name under your personal account. As such, you can commit
   to it as you please without disrupting other developers.
#. Clone locally your fork of F4Enix and create a new branch off of the ``developing`` one.
#. Make your changes on the new branch that you intend to have included in ``developing``.
#. Issue a pull request from GitHub and select the ``developing`` branch of F4Enix main
   repo as the target.
   At a minimum, you should describe what the changes you’ve made are and why you are
   making them. If the changes are related to an outstanding issue, make sure it is
   cross-referenced for its resolution to be properly tracked.
#. A committer will review your pull request based on the criteria above. Any issues with
   the pull request can be discussed directly on the pull request page itself.
#. After the pull request has been thoroughly vetted, it is merged back into the develop
   branch of F4Enix main repo.

.. tip::
    When working with Python API it is strongly recommended to install F4Enix in editable
    mode using pip:

    ```
    pip install -e .
    ```

    in this way the package can be imported from a Python interpreter and any changes made
    are immediately reflected in the installed version
    (that is, you don’t need to keep reinstalling it).

.. tip::
    When working with editable package, often code editor will not be able to perform
    correctly autocompletion and docstring preview. This can be solved in general adding
    the path to your local F4Enix repository to the list of paths where the code editor
    looks for the python interpreter. This fix will be editor dependent but one should 
    be able to find enough information on the web on how to solve it.

Requirements for a successful merge
===================================
The following are minimum requirements necessary for the approval of a pull request:

* the python code should adhere to the `PEP 8 <https://peps.python.org/pep-0008/>`_ convention.
  This can be achieved for instance using `pycodestyle <https://pypi.org/project/pycodestyle/>`_
  as linter in your code editor of choice.
* if a new feature is developed, new test cases must be added to unit test suites.
  `pytest <https://docs.pytest.org/en/7.4.x/>`_ must be used.
* no conflicts are allowed with the ``developing`` branch, i.e., the original ``developing`` branch
  should be pulled into the fork and all eventual conflicts resolved prior to the submission
  of the pull request.
* the new code shall not break any pre-existing feature, i.e., all unit tests and regression tests
  are passed.
* if a new feature is added, it should be properly reported in the sphinx documentation.
  According to the entity of the modifications these may include:

  - docstring (numpy style) including examples and attributes for the main classes;
  - adjourn main body of documentation;
  - add examples on the compiled jupyter notebooks.

Modify documentation using Sphinx
=================================

This documentation is written with
`Sphynx <https://www.sphinx-doc.org/en/master/index.html>`_ using a template
provided by `Read The Docs <https://readthedocs.org/>`_. Before attempting
to modify the documentation, the developer should familiarize with these tools
and with the RST language that is used to actually write the doc. 

Inside the ``docs`` folder of F4Enix repo are located the *source* and *build* directories
of the documentation. To apply a modification, the user must simply modify/add one
or more files in the *source* tree and in the *docs* folder execute from terminal
the ``make html`` command to check that compilation works as intended.

Even if the documentation is not rebuilt locally, a new version is automatically
compiled by ReadTheDocs every time is performed a push to the main branch 
(similarly to what happens with automatic testing of the code).
