.. _installation:

Installation
============

Dipper requires Python version 3.5 or higher

Install with pip:

::

    pip install dipper

Development version
-------------------

The development version can be pulled from GitHub.

::

    pip3 install git+git://github.com/monarch-initiative/dipper.git


Building locally
-------------------

::

    git clone https://github.com/monarch-initiative/dipper.git
    cd dipper
    pip install .

Alternatively, a subset of source specific requirements may be downloaded.  To download the core requirements:

::

    pip install -r requirements.txt

To download source specific requirements use the requirements/ directory, for example:

::

    pip install -r requirements/mgi.txt

To download requirements for all sources:

::

    pip install -r requirements/all-sources.txt