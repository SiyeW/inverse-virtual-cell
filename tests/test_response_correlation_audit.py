import hashlib
import io
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_rna_response_correlations import BudgetReader, calculate, corr, matched_other, run, prefix_array, run_example


def fixture():
    reference = np.array([1., 2, 4, 8, 16, 32])
    delta1 = np.array([.1, -.1, .1, -.1, .1, -.1])
    delta2 = np.array([.1, .1, -.1, -.1, .1, .1])
    means = np.array([reference, reference, (reference+1)*2**delta1-1,
                      (reference+1)*2**delta1-1, (reference+1)*2**delta2-1])
    records = pd.DataFrame(dict(guide=['n1','n2','a1','a2','b1'], target=['NTC','NTC','A','A','B'],
        source_control_type=['non_targeting_control']*2+['targeting']*3,
        cells=[100,100,40,30,30], direct_feature_index=[np.nan]*5))
    effects = np.log2((means+1)/(reference+1))
    pairs = pd.DataFrame([dict(target='A',guide_a='a1',guide_b='a2',cells_a=40,cells_b=30,pearson=1.)])
    return records,pairs,means,effects,np.array(records.guide,dtype=str),np.array(list('abcdef'))


def test_distinguishes_baseline_and_response_and_matches_coverage():
    args = fixture()
    table,example,summary = calculate(*args)
    full = table.iloc[0]
    assert full.response_r == pytest.approx(1)
    assert full.different_target_guide == 'b1' and full.different_target_cells == 30
    assert abs(full.different_target_r) < .01
    assert summary['example']['guide_a']=='a1' and len(example)==6
    assert summary['by_filter']['ntc_mean_ge_10']['features']==2


def test_high_baseline_similarity_does_not_imply_matching_changes():
    _,_,means,effects,_,_ = fixture()
    assert corr(means[2],means[4]) > .99
    assert abs(corr(effects[2],effects[4])) < .01
    assert corr(np.ones(6),np.arange(6)) is None


def test_rejects_row_order_or_effect_miscalculation():
    args = list(fixture())
    args[4] = args[4][::-1]
    with pytest.raises(ValueError,match='Guide order'):
        calculate(*args)
    args = list(fixture());args[3]=args[3].copy();args[3][2,0]+=.1
    with pytest.raises(AssertionError):
        calculate(*args)


def test_direct_feature_mask_must_match():
    args=list(fixture());args[0].loc[2,'direct_feature_index']=1
    with pytest.raises(ValueError,match='mask differs'):
        calculate(*args)
    args[3][2,1]=np.nan
    calculate(*args)


def test_io_budget_and_deadline_fail_without_unbounded_read():
    reader=BudgetReader(io.BytesIO(b'abcdefgh'),4,time.monotonic()+10)
    assert reader.read(3)==b'abc'
    with pytest.raises(ValueError,match='budget'):
        reader.read(2)
    reader=BudgetReader(io.BytesIO(b'a'),4,time.monotonic()-1)
    with pytest.raises(TimeoutError):reader.read(1)


def test_prefix_reader_consumes_only_requested_numeric_rows():
    import zipfile
    values=np.arange(1200,dtype=float).reshape(100,12)
    payload=io.BytesIO();np.savez_compressed(payload,means=values)
    payload.seek(0)
    with zipfile.ZipFile(payload) as archive:
        prefix,shape=prefix_array(archive,'means',3)
    assert shape==(100,12)
    np.testing.assert_array_equal(prefix,values[:3])


def test_accepted_fixture_end_to_end_no_source_or_output_overwrite(tmp_path):
    records,pairs,means,effects,guides,genes=fixture()
    core=tmp_path/'core';(core/'tables').mkdir(parents=True)
    records.to_csv(core/'tables/guide_effect_metadata.csv',index=False)
    pairs.to_csv(core/'tables/within_target_pairs.csv',index=False)
    np.savez_compressed(core/'tables/rna_effect_profiles.npz',means=means,log2fc=effects,guides=guides,gene_ids=genes)
    sources={p.name:p.read_bytes() for p in (core/'tables').iterdir()}
    manifest={f'tables/{p.name}':dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (core/'tables').iterdir()}
    raw=json.dumps(manifest).encode();(core/'output_manifest.json').write_bytes(raw)
    digest=hashlib.sha256(raw).hexdigest();out=tmp_path/'audit'
    result=run(core,out,digest,'fixture',1,10)
    assert result['stage']=='complete' and result['compressed_profile_bytes_read']<1024**2
    assert sources=={p.name:p.read_bytes() for p in (core/'tables').iterdir()}
    with pytest.raises(ValueError,match='unused'):
        run(core,out,digest,'fixture',1,10)


def test_example_only_end_to_end_reports_partial_scope_and_crc_limit(tmp_path):
    records,pairs,means,effects,guides,genes=fixture()
    records.loc[records.target.eq('A'),'source_control_type']='positive_control'
    order=[2,3,4,0,1];records=records.iloc[order].reset_index(drop=True)
    means,effects,guides=means[order],effects[order],guides[order]
    core=tmp_path/'core';(core/'tables').mkdir(parents=True)
    records.to_csv(core/'tables/guide_effect_metadata.csv',index=False)
    pairs.to_csv(core/'tables/within_target_pairs.csv',index=False)
    np.savez_compressed(core/'tables/rna_effect_profiles.npz',means=means,log2fc=effects,guides=guides,gene_ids=genes,
                        reference_means=np.tile(np.array([1.,2,4,8,16,32]),(len(records),1)))
    manifest={f'tables/{p.name}':dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (core/'tables').iterdir()}
    raw=json.dumps(manifest).encode();(core/'output_manifest.json').write_bytes(raw)
    result=run_example(core,tmp_path/'example',hashlib.sha256(raw).hexdigest(),'fixture',1,10)
    assert result['target']=='A' and result['prefix_rows']==2 and result['pairs_reproduced']==1
    assert not result['partial_numeric_members_crc_checked']
    assert not result['full_different_target_background_complete']
    assert result['by_filter']['all_accepted']['response']['median']==pytest.approx(1)
