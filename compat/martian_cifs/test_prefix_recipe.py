import importlib.util,sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).parent))
spec=importlib.util.spec_from_file_location('prefix_replay',Path(__file__).with_name('replay_prefix_fixture.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

@pytest.mark.parametrize('bad',['none','full_run','too_many_pairs','raw_path','missing_crispr'])
def test_prefix_only_guard(bad):
    status={'stage':'real_library_prefix_engineering_passed','scientific_result':False,'forced_cells':10}
    libraries=[{'fastq_dir':'/prefix/rna','fastq_id':'rna','feature_type':'Gene Expression'},{'fastq_dir':'/prefix/crispr','fastq_id':'crispr','feature_type':'CRISPR Guide Capture'}]
    receipt={'scientific_result':False,'pairs_per_library':100000,'libraries':libraries}
    rows=[{'fastqs':x['fastq_dir'],'sample':x['fastq_id'],'library_type':x['feature_type']} for x in libraries]
    if bad=='full_run':status['scientific_result']=True
    if bad=='too_many_pairs':receipt['pairs_per_library']=100000000
    if bad=='raw_path':rows[0]['fastqs']='/raw/full'
    if bad=='missing_crispr':receipt['libraries']=libraries[:1]
    if bad=='none':module.validate_recipe(status,receipt,rows)
    else:
        with pytest.raises(ValueError):module.validate_recipe(status,receipt,rows)
