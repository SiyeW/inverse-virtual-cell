import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import pytest
from compare_energy_fixed_n import fixed_rows,same_cohort_baseline


def test_fixed_draw_invariant_to_row_order_and_preserves_members():
    bc=np.array(['a','b','c','d','e','f']);rows=np.arange(6)
    a=fixed_rows(rows,bc,'g',3,246);b=fixed_rows(rows[::-1],bc,'g',3,246)
    assert np.array_equal(a,b) and len(set(a))==3 and set(a)<=set(rows)
    assert np.array_equal(fixed_rows(rows,bc,'g',6,246),fixed_rows(rows,bc,'g',6,246))
    with pytest.raises(ValueError):fixed_rows(rows,bc,'g',7,246)
    with pytest.raises(ValueError):fixed_rows([0,0],bc,'g',2,246)


def test_baseline_correction_recomputed_only_on_same_cohort():
    f=pd.DataFrame(dict(guide=['a','b','c'],pvalue=[.01,.04,.9],pvalue_adj=[.03,.06,.9]))
    out=same_cohort_baseline(f,['b','a'])
    assert out.guide.tolist()==['b','a']
    assert np.allclose(out.pvalue_adj,[.04,.02]) and out.significant_adj.all()
    with pytest.raises(ValueError):same_cohort_baseline(f,['a','a'])
