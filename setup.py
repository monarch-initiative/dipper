#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(name='Dipper',
      version='0.0.1',
      description='Data Ingest Pipeline',
      packages=find_packages(),
      install_requires=['psycopg2', 'rdflib', 'isodate', 'roman', 'python-docx', 'pyyaml', 'pysftp', 'biopython',
                        'docx', 'beautifulsoup4'],
      include_package_data=True
      )
