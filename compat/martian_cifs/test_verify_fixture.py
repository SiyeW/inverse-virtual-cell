import importlib.util
from pathlib import Path
import h5py
import pytest
import json

spec=importlib.util.spec_from_file_location('verify_fixture',Path(__file__).with_name('verify_fixture.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

@pytest.mark.parametrize('difference',['none','count','barcode','coverage'])
def test_matrix_equivalence(tmp_path,difference):
    a,b=tmp_path/'a.h5',tmp_path/'b.h5'
    for p in [a,b]:
        with h5py.File(p,'w') as f:
            f['matrix/data']=[1,2 if p==a or difference!='count' else 3]
            f['matrix/barcodes']=[b'A',b'C' if p==b and difference=='barcode' else b'B']
            if p==a or difference!='coverage':f['matrix/shape']=[2,2]
    if difference=='none':assert len(module.compare_matrix(a,b))==3
    else:
        with pytest.raises(ValueError):module.compare_matrix(a,b)

@pytest.mark.parametrize('failure',['publication','json','status'])
def test_fixture_rejects_incomplete_acceptance(tmp_path,failure):
    root=tmp_path/'fixture';root.mkdir()
    (root/'status.json').write_text(json.dumps({'stage':'failed_preserved' if failure=='status' else 'engineering_fixture_passed','scientific_result':False}))
    run=root/'tiny_no_bam';(run/'SC_RNA_COUNTER_CS/fork0').mkdir(parents=True)
    (run/'_log').write_text('Could not move output files' if failure=='publication' else '')
    (run/'SC_RNA_COUNTER_CS/fork0/_outs').write_text('{invalid' if failure=='json' else '{}')
    with pytest.raises(ValueError):module.verify(root,root)
