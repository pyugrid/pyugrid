import os
import os.path


def get_test_file_path(file_path):
    """translates a file path to be relative to the test files directory"""
    return os.path.join(os.path.dirname(__file__), 'files', file_path)