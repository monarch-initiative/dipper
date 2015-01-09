__author__ = 'condit@sdsc.edu'


class CurieUtil:

    # curie_map should be in the form URI -> prefix:
    # ie: 'http://foo.org/bar_' -> 'bar'
    # this class will not currently handle multiple uris with different prefixes
    def __init__(self, curie_map):
        self.curie_map = curie_map
        self.uri_map = {}
        for key, value in curie_map.items():
            self.uri_map[value] = key
        return

    # Get a CURIE from a uri
    def get_curie(self, uri):
        for key, value in self.uri_map.items():
            if uri.startswith(key):
                return '%s:%s' % (value, uri[len(key):len(uri)])
        return None

    # Get a URI from a CURIE
    def get_uri(self, curie):
        parts = curie.split(':')
        if 1 == len(parts):
            if (curie != ''):
                print ("ERROR: Not a properly formed curie: \"",curie,"\"",sep='')
            return None
        prefix = parts[0]
        if prefix in self.curie_map:
            return '%s%s' % (self.curie_map.get(prefix), curie[(curie.index(':') + 1):])
        print("ERROR: Curie prefix not defined for",curie)
        return None
