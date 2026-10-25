from pathlib import Path
import sys
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_ntc_reference import calculate,run,shift_stats
from audit_rna_response_correlations import corr


def inputs():
    reference=np.array([1.,2,4,8,16,32])
    ntc=np.array([reference*.8,reference*.9,reference*1.1,reference*1.2])
    counts=np.array([40,30,10,20,30,40])
    center=np.average(ntc,axis=0,weights=counts[2:])
    delta=np.array([.1,-.1,.2,-.2,.1,-.1])
    mean=(center+1)*2**delta-1
    means=np.vstack([mean,mean,ntc])
    records=pd.DataFrame(dict(guide=['a1','a2','n1','n2','n3','n4'],target=['A','A','NTC','NTC','NTC','NTC'],
        source_control_type=['positive_control']*2+['non_targeting_control']*4,cells=counts,
        direct_feature_index=[np.nan]*6))
    pairs=pd.DataFrame([dict(target='A',guide_a='a1',guide_b='a2',cells_a=40,cells_b=30,pearson=1.)])
    return records,means,center,pairs


def test_reference_weighting_leaveout_and_complete_disjoint_partitions():
    records,means,center,pairs=inputs()
    guide,refs,changes,genes,plans,summary=calculate(records,means,center,pairs,len(records),splits=4)
    assert summary['ntc_cells']==100 and summary['ntc_guides']==4
    assert summary['largest_guide_cell_weight']==.4
    assert guide.reference_overlap_cells.eq(0).all()
    first=guide.loc[guide.guide.eq('n1')].iloc[0]
    expected=np.average(means[3:],axis=0,weights=records.cells.to_numpy()[3:])
    assert first.reference_cells==90
    assert first.leave_one_rms_log2_shift==pytest.approx(shift_stats(expected,center)['rms_log2_shift'])
    for split in plans:
        assert not set(split['ntc_guides_a'])&set(split['ntc_guides_b'])
        assert set(split['ntc_guides_a'])|set(split['ntc_guides_b'])==set(records.guide[2:])
        assert split['cells_a']+split['cells_b']==100
    assert not summary['true_biological_replicates'] and not summary['baseline_replacement_adopted']


def test_identical_ntc_guides_give_no_reference_sensitivity():
    records,means,_,pairs=inputs()
    means[2:]=np.array([1.,2,4,8,16,32])
    center=means[2].copy()
    delta=np.array([.1,-.1,.2,-.2,.1,-.1]);means[:2]=(center+1)*2**delta-1
    guide,refs,changes,genes,plans,summary=calculate(records,means,center,pairs,len(records))
    np.testing.assert_allclose(refs.rms_log2_shift,0,atol=1e-14)
    np.testing.assert_allclose(changes.response_r,1,atol=1e-12)


def test_mismatching_pooled_reference_or_old_pair_fails():
    records,means,center,pairs=inputs()
    with pytest.raises(AssertionError):calculate(records,means,center+.1,pairs,len(records))
    pairs.loc[0,'pearson']=.5
    with pytest.raises(ValueError,match='pair differs'):calculate(records,means,center,pairs,len(records))


def test_common_reference_error_can_add_correlation_instead_of_destroy_it():
    # Orthogonal centered error vectors: no simulated data or asymptotic assumption needed.
    a=np.array([1,1,1,1,-1,-1,-1,-1.])
    b=np.array([1,1,-1,-1,1,1,-1,-1.])
    c=np.array([1,-1,1,-1,1,-1,1,-1.])
    d=np.array([1,-1,-1,1,-1,1,1,-1.])
    assert corr(a-c,b-c)==pytest.approx(.5)
    assert corr(a-c,b-d)==pytest.approx(0)


def test_tiny_accepted_source_end_to_end_stays_read_only(tmp_path):
    records,means,center,pairs=inputs();core=tmp_path/'core';(core/'tables').mkdir(parents=True)
    records.to_csv(core/'tables/guide_effect_metadata.csv',index=False)
    pairs.to_csv(core/'tables/within_target_pairs.csv',index=False)
    np.savez_compressed(core/'tables/rna_effect_profiles.npz',means=means,reference_means=np.tile(center,(len(records),1)),
                        guides=records.guide.to_numpy(dtype=str),gene_ids=np.array(list('abcdef')))
    saved={p.name:p.read_bytes() for p in (core/'tables').iterdir()}
    manifest={f'tables/{p.name}':dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (core/'tables').iterdir()}
    raw=json.dumps(manifest).encode();(core/'output_manifest.json').write_bytes(raw)
    result=run(core,tmp_path/'audit',hashlib.sha256(raw).hexdigest(),'fixture',1,10)
    assert result['stage']=='complete' and result['covered_response_pairs']==1
    assert not result['partial_numeric_member_crc_checked']
    assert saved=={p.name:p.read_bytes() for p in (core/'tables').iterdir()}
