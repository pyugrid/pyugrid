"""
assorted utilities useful for the tests
"""

import os
import contextlib


@contextlib.contextmanager
def chdir(dirname=None):
    curdir = os.getcwd()
    try:
        if dirname is not None:
            os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)

def get_test_file_path(file_path):
    """translates a file path to be relative to the test files directory"""
    return os.path.join(os.path.dirname(__file__), 'files', file_path)