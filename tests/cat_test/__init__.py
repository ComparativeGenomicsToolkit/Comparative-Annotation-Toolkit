"""support functions for CAT tests using pytest"""

import os
import os.path as osp
import pytest
import difflib

def get_test_id():
    """get id for the current test being run"""
    return os.environ.get('PYTEST_CURRENT_TEST')

def get_test_cid():
    """get cleaned up id for current test suitable for creating output files"""
    bn = osp.basename(get_test_id())
    return bn.replace(" (call)", "")

def get_test_dir():
    """return directory containing current to find input and output files"""
    # FIXME: change to find the directory to allow running tests
    # from high level directory
    return "."

def get_input_file(fname):
    """return path to file in the test input directory"""
    return osp.join(get_test_dir(), "input", fname)

def get_output_dir():
    """get the path to the output directory to use for this test, create if it doesn't exist"""
    d = osp.join(get_test_dir(), "output")
    os.makedirs(d, exist_ok=True)
    return d

def get_output_file(ext):
    """Get path to the output file, using the current test id and append
    ext, which should contain a dot."""
    f = osp.join(get_output_dir(), get_test_cid() + ext)
    return f

def get_expected_file(ext, basename=None):
    """Get path to the expected file, using the current test id and append
    ext. If basename is used, it is instead of the test id, allowing share
    an expected file between multiple tests."""
    if basename is None:
        basename = get_test_cid()
    return osp.join(get_test_dir(), "expected", basename + ext)

def _get_lines(file):
    with open(file) as fh:
        return fh.readlines()

def diff_files(exp_file, out_file):
    """diff expected and output files."""

    exp_lines = _get_lines(exp_file)
    out_lines = _get_lines(out_file)

    diff = difflib.unified_diff(exp_lines, out_lines, exp_file, out_file)
    cnt = 0  # count because diff is a generator
    for l in diff:
        print(l, end=' ')
        cnt += 1
    assert cnt == 0, f"{exp_file} and {out_file} differ"

def diff_expected(ext):
    """diff expected and output files, with names computed from test id."""
    diff_files(get_expected_file(ext), get_output_file(ext))
