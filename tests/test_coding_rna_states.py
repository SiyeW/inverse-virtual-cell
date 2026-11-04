from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).parents[1]/"scripts"))
from analyze_coding_rna_states import comparison_masks,residualize_amount,profile_correlation,all_profile_pairs


def positive():
    return dict(frameshift_max_same_variant_families=1,indel_families=1)


def reference():
    return dict(frameshift_site_reference_max=1)


def test_pure_guide_records_are_not_compared_against_another_guides_evidence():
    rows=[positive(),positive(),reference(),reference(),positive(),reference(),{}]
    guides=np.array(['a','a','b','b','c','c','c'])
    pos,neg,use,paired=comparison_masks(rows,guides,'internal_frameshift')
    assert paired==['c']
    assert np.array_equal(use,[0,0,0,0,1,1,0])
    assert pos.sum()==3 and neg.sum()==3


def test_empty_evidence_does_not_fabricate_a_comparison():
    p,n,use,paired=comparison_masks([],np.array([],dtype=str),'internal_frameshift')
    assert len(p)==len(n)==len(use)==0 and paired==[]


def test_RNA_amount_adjustment_does_not_use_labels_and_removes_linear_amount_term():
    total=np.array([9,19,29,39,49,59]);z=np.log1p(total)
    x=np.column_stack([2+3*z,5-z,z*z])
    residual=residualize_amount(x,total)
    assert np.max(abs(residual[:,:2]))<1e-12
    assert np.max(abs(np.column_stack([np.ones(6),z]).T@residual))<1e-11


def test_pairwise_RNA_universe_excludes_both_directly_targeted_RNAs():
    genes=np.array(['targetA','targetB','x','y','z'])
    left=np.array([1000,-1000,1,2,3]);right=np.array([-900,800,4,5,6])
    assert np.isclose(profile_correlation(left,right,genes,'targetA','targetB'),1)


def test_undefined_flat_profiles_are_not_zero_correlations():
    assert np.isnan(profile_correlation(np.ones(3),np.arange(3),np.array(['x','y','z']),'a','b'))


def test_every_same_target_and_other_target_pair_is_used_without_sampling():
    profiles=[]
    for i,target in enumerate(['a','a','a','b','b']):
        profiles.append(dict(guide=str(i),target_gene_id=target,n_positive=1+i,n_reference=3,
                             delta=np.array([0.,float(i+1),2.]),delta_amount_adjusted=np.array([0.,float(i+1),2.])))
    pairs=all_profile_pairs(profiles,np.array(['x','y','z']),'company')
    assert len(pairs)==10
    assert sum(r['same_target'] for r in pairs)==4  # C(3,2)+C(2,2)
    assert len({(r['guide_left'],r['guide_right']) for r in pairs})==10
