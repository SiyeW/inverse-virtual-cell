from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from diagnose_control_guides import centroid_null_test


def test_null_constant_and_separated_groups():
    effect, p = centroid_null_test(np.zeros((20,3)),np.zeros((80,3)),draws=99)
    assert effect == 0 and p == 1
    effect, p = centroid_null_test(np.ones((20,3))*10,np.zeros((80,3)),draws=99)
    assert effect > 10 and p == .01
