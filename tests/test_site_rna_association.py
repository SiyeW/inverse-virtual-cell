import sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from check_site_rna_association import energy_from_indicators, within_guide_permutations, adjusted_bh


def test_energy_matches_direct_group_distance_formula():
    x = np.array([[0., 0], [2, 1], [1, 0], [5, 2], [3, 2]])
    labels = np.array([1, 1, 0, 0, 0]); d = cdist(x, x)
    a, b = labels == 1, labels == 0
    expected = 2*d[np.ix_(a,b)].mean()-d[np.ix_(a,a)].sum()/2-d[np.ix_(b,b)].sum()/6
    assert np.allclose(energy_from_indicators(d, labels), expected)


def test_permutations_preserve_counts_within_each_sgrna():
    labels = np.array([1,1,0,0,1,0,0]); guide = np.array(['a']*4+['b']*3)
    p, mutable = within_guide_permutations(labels, guide, repeats=50)
    assert mutable == 7
    assert np.all(p[:,:4].sum(axis=1)==2) and np.all(p[:,4:].sum(axis=1)==1)


def test_no_cross_guide_shuffle_when_labels_are_fixed_per_guide():
    y = np.array([1,1,0,0]); g = np.array(['a','a','b','b'])
    p, mutable = within_guide_permutations(y,g,repeats=10)
    assert mutable == 0 and np.array_equal(p,np.repeat(y[None,:],10,axis=0))


def test_bh_retains_order_and_known_values():
    assert np.allclose(adjusted_bh([.04,.001,.03]), [.04,.003,.04])
