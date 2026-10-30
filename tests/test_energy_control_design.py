from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pytest
from check_energy_control_design import random_query
from diagnose_energy_reference import independent_reference


def test_null_uses_control_population_and_preserves_disjoint_reference():
    labels=np.array(['own','own','a','a','b','b','target']);roles=np.array(['non_targeting_control']*6+['targeting']);bc=np.array(list('abcdefg'))
    ref=np.array([2,4]);q=random_query('own',labels,roles,bc,ref,2,812)
    assert set(q)=={3,5} and not set(q).intersection(ref)
    assert np.array_equal(q,random_query('own',labels,roles,bc,ref,2,812))
    with pytest.raises(ValueError):random_query('own',labels,roles,bc,ref,3,812)


def test_independent_reference_excludes_entire_own_guide():
    labels=np.array(['own','own','a','a','b','b','target']);roles=np.array(['non_targeting_control']*6+['targeting']);bc=np.array(list('abcdefg'))
    ref=independent_reference('own',labels,roles,bc,3,123)
    assert len(set(ref))==3 and set(ref)<=set([2,3,4,5])
