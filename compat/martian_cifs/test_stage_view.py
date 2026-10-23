import importlib.util
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('stage_view',Path(__file__).with_name('stage_view.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def fixture(tmp_path):
    src=tmp_path/'source';(src/'bin').mkdir(parents=True);(src/'external/martian/bin').mkdir(parents=True)
    (src/'bin/cellranger').write_bytes(b'original scientific launcher')
    (src/'external/martian/bin/mrp').write_bytes(b'original runtime')
    (src/'external/martian/bin/mrjob').write_bytes(b'original adapter')
    try:(src/'cellranger').symlink_to('bin/cellranger')
    except OSError:pytest.skip('Native symlink permission required for software view test')
    runtime=tmp_path/'candidate';runtime.write_bytes(b'candidate runtime')
    return src,runtime

def test_original_unchanged_and_copies_independent(tmp_path):
    src,runtime=fixture(tmp_path);out=tmp_path/'view'
    module.stage(src,runtime,out)
    assert (out/'cellranger').resolve()==out/'bin/cellranger'
    assert (out/'external/martian/bin/mrp').read_bytes()==b'candidate runtime'
    (out/'bin/cellranger').write_bytes(b'test view changed')
    assert (src/'bin/cellranger').read_bytes()==b'original scientific launcher'
    assert (src/'external/martian/bin/mrp').read_bytes()==b'original runtime'

def test_existing_target_refused(tmp_path):
    src,runtime=fixture(tmp_path);out=tmp_path/'view';out.mkdir()
    with pytest.raises(ValueError):module.stage(src,runtime,out)

def test_nested_source_refused(tmp_path):
    src,runtime=fixture(tmp_path)
    with pytest.raises(ValueError):module.stage(src,runtime,src/'view')
