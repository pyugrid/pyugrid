"""
assorted utilities useful for the tests
"""

from __future__ import (absolute_import, division, print_function)

import os
import contextlib


@contextlib.contextmanager
def chdir(dirname=None):
    curdir = os.getcwd()
    try:
        if dirname is not None:
            if dirname == 'files':
                import test
                dirname = os.path.split(test.__file__)[0]
                dirname = os.path.join(dirname, 'files')
            os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)

def get_test_file_path(file_path):
    """translates a file path to be relative to the test files directory"""
    return os.path.join(os.path.dirname(__file__), 'files', file_path)

