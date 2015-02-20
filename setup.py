#!/usr/bin/env python3

from distutils.core import setup

# TODO Move packages into single directory (named dipper?) to
# namespace
setup(name='DIPper',
      version='1.0',
      description='Data Ingest Pipeline',
      packages=['sources', 'models', 'utils'],
      install_requires=['psycopg2', 'rdflib', 'isodate', 'roman', 'python-docx', 'pyyaml']
      )
