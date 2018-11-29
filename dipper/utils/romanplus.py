"""
Convert to and from Roman numerals
This program is part of "Dive Into Python", a free Python tutorial for
experienced programmers.  Visit http://diveintopython.org/ for the
latest version.

This program is free software; you can redistribute it and/or modify
it under the terms of the Python 2.1.1 license, available at
http://www.python.org/2.1.1/license.html

Note:
This has been modified to add optional characters
after the initial roman numbers by nlw.
"""
import re

__author__ = "Mark Pilgrim (f8dy@diveintopython.org)"
__version__ = "1.4"
__date__ = "8 August 2001"
__copyright__ = "Copyright (c) 2001 Mark Pilgrim"

# Define digit's regular expression mapping
romanNumeralMap = (('M', 1000),
                   ('CM', 900),
                   ('D', 500),
                   ('CD', 400),
                   ('C', 100),
                   ('XC', 90),
                   ('L', 50),
                   ('XL', 40),
                   ('X', 10),
                   ('IX', 9),
                   ('V', 5),
                   ('IV', 4),
                   ('I', 1))


def toRoman(num):
    """convert integer to Roman numeral"""
    if not 0 < num < 5000:
        raise ValueError("number %n out of range (must be 1..4999)", num)
    if int(num) != num:
        raise TypeError("decimals %n can not be converted", num)

    result = ""
    for numeral, integer in romanNumeralMap:
        while num >= integer:
            result += numeral
            num -= integer
    return result


# Define pattern to detect valid Roman numerals
romanNumeralPattern = re.compile("""
    ^                   # beginning of string
    M{0,4}              # thousands - 0 to 4 M's
    (CM|CD|D?C{0,3})    # hundreds - 900 (CM), 400 (CD), 0-300 (0 to 3 C's),
                        #            or 500-800 (D, followed by 0 to 3 C's)
    (XC|XL|L?X{0,3})    # tens - 90 (XC), 40 (XL), 0-30 (0 to 3 X's),
                        #        or 50-80 (L, followed by 0 to 3 X's)
    (IX|IV|V?I{0,3})    # ones - 9 (IX), 4 (IV), 0-3 (0 to 3 I's),
                        #        or 5-8 (V, followed by 0 to 3 I's)
    [A-Z]               # optional suffix letter,
                        #        but don't retain.
                        #         differs from original roman.py
    $                   # end of string
    """, re.VERBOSE)


def fromRoman(strng):
    """convert Roman numeral to integer"""
    if not strng:
        raise TypeError('Input can not be blank')
    if not romanNumeralPattern.search(strng):
        raise ValueError('Invalid Roman numeral: %s', strng)

    result = 0
    index = 0
    for numeral, integer in romanNumeralMap:
        while strng[index:index+len(numeral)] == numeral:
            result += integer
            index += len(numeral)
    return result
