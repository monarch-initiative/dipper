from utils.TestUtils import TestUtils
import logging

logger = logging.getLogger(__name__)

def compare_local_remote_bytes(source):
    """
    test to see if fetched file is the same size as the remote file
    using information in the content-length field in the HTTP header
    :return: True or False
    """
    is_equal = True
    test = TestUtils()
    source.rawdir = "../raw/"+source.name
    source.fetch(False)
    for file, paths in source.files.items():
        file_path = '/'.join((source.rawdir, paths['file']))
        local_size = test.get_local_file_size(file_path)
        remote_size = test.get_remote_content_len(paths['url'])
        if local_size != int(remote_size):
            is_equal = False
            logger.warn('local file and remote file different sizes\n'
                        '%s has size %s, %s has size %s', file_path,
                        local_size, paths['url'], remote_size)

    return is_equal
