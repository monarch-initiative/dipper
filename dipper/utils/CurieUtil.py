
import logging

__author__ = 'condit@sdsc.edu'

logger = logging.getLogger(__name__)


class CurieUtil(object):
    '''
    Create compact URI
    '''
    def __init__(self, curie_map):
        '''
        curie_map format is: curie_prefix -> URI_prefix:
        ie: 'bar': 'http://foo.org/bar_'

        '''
        self.curie_map = curie_map
        if curie_map is not None:  # inverse the map
            if len(set(curie_map.keys())) < len(set(curie_map.values())):
                logger.warrning("Curie map is NOT one to one!")
                logger.warrning(
                "`get_curie_prefix(IRI)` will return the same prefix for different base IRI")
            self.uri_map = {}
            for key, value in curie_map.items():
                self.uri_map[value] = key
        return

    def get_curie(self, uri):
        '''Get a CURIE from a URI '''
        prefix = self.get_curie_prefix(uri)
        if prefix is not None:
            key = self.curie_map[prefix]
            return '%s:%s' % (prefix, uri[len(key):len(uri)])
        return None

    def get_curie_prefix(self, uri):
        ''' Return the CURIE's prefix:'''
        for key, value in self.uri_map.items():
            if uri.startswith(key):
                return value
        return None

    def get_uri(self, curie):
        ''' Get a URI from a CURIE '''
        if curie is None:
            return None
        parts = curie.split(':')
        if len(parts) == 1:
            if curie != '':
                logger.error("Not a properly formed curie: \"%s\"", curie)
            return None
        prefix = parts[0]
        if prefix in self.curie_map:
            return '%s%s' % (self.curie_map.get(prefix),
                             curie[(curie.index(':') + 1):])
        logger.error("Curie prefix not defined for %s", curie)
        return None

    def prefix_exists(self, pfx):
        return pfx in self.curie_map
