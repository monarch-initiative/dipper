.. _configuration:

Configuring dipper with keys and passwords
==========================================

Add private configuration parameters into your private conf.yaml file.
Examples of items to put into the config include:

* database connection parameters (in the "dbauth" object)
* ftp login credentials
* api keys (in the "keys" object)

These are organized such that within any object (dbauth, keys, etc),
they are keyed again by the source's name.

Here is an example::
{
   "keys": {
      "omim": "myomimkey1234",
      "omimftp": "myomimftpkey5678",
      "ncbi": "myncbikey",
   },
   "dbauth": {
      "mgi": {
         'user': 'mymgiuser',
         'password': 'mymgipw',
         'host': 'mgi-adhoc.jax.org',
         'database': 'mgd',
         'port': 5432
      }
   }
}

This file must be placed in the dipper package directory and named conf.yaml.
If building locally this is in the dipper/dipper/ directory. If installed with pip
this will be in path/to/env/lib/python3.x/site-packages/dipper/ directory.
