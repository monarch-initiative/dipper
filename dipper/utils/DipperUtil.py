__author__ = 'nlw'


import logging
import re
import unicodedata

logger = logging.getLogger(__name__)


class DipperUtil:
    """
    Various utilities and quick methods used in this application
    """

    def flatten(self, l):
        """
        Remove None from an array or list
        :param l: An array
        :return:  An array with the None elements removed
        """
        l = list(filter(None.__ne__, l))

        return l

    def remove_control_characters(self, s):
        return "".join(ch for ch in s if unicodedata.category(ch)[0] != "C")


