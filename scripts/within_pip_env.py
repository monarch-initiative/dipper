'''
    Determine if a shell is within a python virtual enviroment or not
    mimics program return code "$?" on stdout without the stderr cruft on failure
    which simplifies assignment within enviroments such as make

    hence  '0' is no error (detecting pip venv exists)
    and    '1' means pip venv is not detectable due to some error (e.g. not existing)

    '''

import sys

x = '0'
try:
    sys.real_prefix
except:
    x = '1'

sys.stdout.write(x)
sys.stdout.flush()
