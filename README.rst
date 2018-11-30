===================
HLA in the Pipeline
===================


.. image:: https://img.shields.io/pypi/v/hlapipeline.svg
        :target: https://pypi.python.org/pypi/hlapipeline

.. image:: https://img.shields.io/travis/stsci-hack/hlapipeline.svg
        :target: https://travis-ci.org/stsci-hack/hlapipeline

.. image:: https://readthedocs.org/projects/hlapipeline/badge/?version=latest
        :target: https://hlapipeline.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Code for implementing HLA-type processing in the HST Pipeline


* Free software: BSD license
* Documentation: https://hlapipeline.readthedocs.io.


Installation
------------

* TODO

Running Tests
--------------
This package comes with a number of tests which can be run using data from the HST
archive as accessed using `astroquery`.  These tests can be run after installing
this package using the following steps:

  * Define the environment variable `TEST_BIGDATA` which is the location of input
    test data.  It can either be:

      * a local directory (such as `/grp/hst/ssb/artifactory/`) or
      * an Artifactory server (such as `https://bytesalad.stsci.edu/artifactory`)

  * (Optional) change to the tests directory in the hlapipeline package source code
  * Run pytest on the desired tests.  For example::

    pytest -s --basetemp=/internal/1/pytest-hst --bigdata test_align.py >& test_align.log

    where the parameters are defined as:

      * `-s`: disable all capturing of stdout/stderr within pytest so that it will
        will show up interactively
      * `--basetemp=`: this defines what LOCAL directory will be used for running
        the tests and writing out the results.  All contents of this directory
        will be deleted upon the next run of the tests.

.. warning ::
  Several of these tests require the use of a large number of datasets and therefore
  not only take a lot of disk space to run, but also can take a significant amount
  of time to run.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
