'''
    Replace me
    
'''

class Template(Source):
    '''
        Replace me
    '''
    files = {
        'src_key_1': {
            'file': 'filename.ext',
            'url':  'protocal://full/path/filename.ext',
            'headers': {'User-Agent': USER_AGENT},
            'clean': 'some_url',  # what to say when you do not want to say much
            'columns': (  # expected
                'first_column_name',
                'second_column_name',
                'third_column_name',
            ),
    }

    
    def __init__(
        self,
        graph_type='rdf_graph',     # or streamed_graph
        are_bnodes_skized=True,     # typically True
        name=None,                  # identifier (lowercase ingest file name
        ingest_title=None,
        ingest_url=None,
        license_url=None,           # _only_ if it is _our_ lic
        data_rights=None,           # external page that points to their current lic
        file_handle=None
    ):
        '''
            Replace me
        '''
        # super() / parent Source class  resources include:
        #   curie_map       dict
        #   all_test_ids    dict
        #   globaltt        dict  and its inverse globaltcid
        #   localtt         dict
        #   resolve()       function over local & global tt   
        
        return
    
    # 
    def fetch(self, is_dl_forced=False, files=None):
        '''
            Replace me
        '''
        self.get_files(is_dl_forced, files)
        return

    def parse(self):
        '''
            Replace me
        '''
        
        return

    def process_src_key_1():
        
    
