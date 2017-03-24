#!/usr/bin/env python3

import re

__author__ = 'Mahmoud Adel <mahmoud.adel2@gmail.com>'
__version__ = 0.4

# from https://raw.githubusercontent.com/mahmoudadel2/pysed/master/pysed.py
# with a few edits from us


def replace(oldstr, newstr, infile, dryrun=False):
    """
    Sed-like Replace function..
    Usage: pysed.replace(<Old string>, <Replacement String>, <Text File>)
    Example: pysed.replace('xyz', 'XYZ', '/path/to/file.txt')

    This will dump the output to STDOUT instead of changing the input file.
    Example 'DRYRUN':
    pysed.replace('xyz', 'XYZ', '/path/to/file.txt', dryrun=True)

    """

    linelist = []
    with open(infile) as f:
        for item in f:
            newitem = re.sub(oldstr, newstr, item)
            linelist.append(newitem)
    if dryrun is False:
        with open(infile, "w") as f:
            f.truncate()
            for line in linelist:
                f.writelines(line)
    elif dryrun is True:
        for line in linelist:
            print(line, end='')
    else:
        exit("""Unknown option specified to 'dryrun' argument,
             Usage: dryrun=<True|False>.""")


def rmlinematch(oldstr, infile, dryrun=False):
    """
    Sed-like line deletion function based on given string..
    Usage: pysed.rmlinematch(<Unwanted string>, <Text File>)
    Example: pysed.rmlinematch('xyz', '/path/to/file.txt')
    Example:
    'DRYRUN': pysed.rmlinematch('xyz', '/path/to/file.txt', dryrun=True)
    This will dump the output to STDOUT instead of changing the input file.

    """

    linelist = []
    with open(infile) as f:
        for item in f:
            rmitem = re.match(r'.*{}'.format(oldstr), item)
            if isinstance(rmitem) == isinstance(None):
                linelist.append(item)
    if dryrun is False:
        with open(infile, "w") as f:
            f.truncate()
            for line in linelist:
                f.writelines(line)
    elif dryrun is True:
        for line in linelist:
            print(line, end='')
    else:
        exit("""Unknown option specified to 'dryrun' argument,
              Usage: dryrun=<True|False>.""")


def rmlinenumber(linenumber, infile, dryrun=False):
    """
    Sed-like line deletion function based on given line number..
    Usage: pysed.rmlinenumber(<Unwanted Line Number>, <Text File>)
    Example: pysed.rmlinenumber(10, '/path/to/file.txt')
    Example 'DRYRUN': pysed.rmlinenumber(10, '/path/to/file.txt', dryrun=True)
    #This will dump the output to STDOUT instead of changing the input file.

    """

    linelist = []
    linecounter = 0
    if isinstance(linenumber) != isinstance(linecounter):
        exit("""'linenumber' argument must be an integer.""")
    with open(infile) as f:
        for item in f:
            linecounter = linecounter + 1
            if linecounter != linenumber:
                linelist.append(item)
    if dryrun is False:
        with open(infile, "w") as f:
            f.truncate()
            for line in linelist:
                f.writelines(line)
    elif dryrun is True:
        for line in linelist:
            print(line, end='')
    else:
        exit("""Unknown option specified to 'dryrun' argument,
              Usage: dryrun=<True|False>.""")
