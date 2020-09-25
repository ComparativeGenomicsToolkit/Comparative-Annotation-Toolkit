"""tests of parent gene assignment"""

import sys
import pytest
sys.path = ["..", "../.."] +  sys.path
from cat_test import get_input_file, get_output_file, diff_expected
from parent_gene_assignment import assign_parents


def test_parent_assign_chm13():
    """parent assignment cases from T2T CHM13 20200727 assembly"""
    pass
