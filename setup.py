#!/usr/bin/env python

from distutils.core import setup

setup(name='DIPper',
      version='1.0',
      description='Data Ingest Pipeline',
      packages=['dipper', 'dipper.sources', 'dipper.models', 'dipper.utils'],
      package_dir={'dipper':'/', 'dipper.sources': '/',
                   'dipper.models':'/', 'dipper.utils':'/'
                   },
      install_requires=['psycopg2', 'rdflib', 'isodate', 'roman', 'python-docx', 'pyyaml']
      )
