import importlib.util
from pathlib import Path
import pytest
import json

spec=importlib.util.spec_from_file_location('stage_view',Path(__file__).with_name('stage_view.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def fixture(tmp_path):
    src=tmp_path/'source';(src/'bin').mkdir(parents=True);(src/'external/martian/bin').mkdir(parents=True)
    (src/'bin/cellranger').write_bytes(b'original scientific launcher')
    (src/'external/martian/bin/mrp').write_bytes(b'original runtime')
    (src/'external/martian/bin/mrjob').write_bytes(b'original adapter')
    (src/'.env.json').write_text(json.dumps({'envs':[{'key':'LD_PRELOAD','type':'unset','value':''},{'key':'KEEP','type':'string','value':'unchanged'}]}))
    try:(src/'cellranger').symlink_to('bin/cellranger')
    except OSError:pytest.skip('Native symlink permission required for software view test')
    runtime=tmp_path/'candidate';runtime.write_bytes(b'candidate runtime')
    shim=tmp_path/'shim';shim.write_bytes(b'candidate shim')
    return src,runtime,shim

def test_original_unchanged_and_copies_independent(tmp_path):
    src,runtime,shim=fixture(tmp_path);out=tmp_path/'view'
    module.stage(src,runtime,out,shim)
    assert (out/'cellranger').resolve()==out/'bin/cellranger'
    assert (out/'external/martian/bin/mrp').read_bytes()==b'candidate runtime'
    (out/'bin/cellranger').write_bytes(b'test view changed')
    assert (src/'bin/cellranger').read_bytes()==b'original scientific launcher'
    assert (src/'external/martian/bin/mrp').read_bytes()==b'original runtime'
    original=json.loads((src/'.env.json').read_text());derived=json.loads((out/'.env.json').read_text())
    assert original['envs'][0]['type']=='unset'
    assert derived['envs'][0]=={'key':'LD_PRELOAD','type':'string','value':str(out/'compat/same_mode_fchmod.so')}
    assert original['envs'][1:]==derived['envs'][1:]

def test_existing_target_refused(tmp_path):
    src,runtime,shim=fixture(tmp_path);out=tmp_path/'view';out.mkdir()
    with pytest.raises(ValueError):module.stage(src,runtime,out,shim)

def test_nested_source_refused(tmp_path):
    src,runtime,shim=fixture(tmp_path)
    with pytest.raises(ValueError):module.stage(src,runtime,src/'view',shim)
