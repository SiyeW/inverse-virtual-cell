import importlib.util
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('inspect_clone',Path(__file__).with_name('clone_metadata_for_inspect.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def source(tmp_path):
    src=tmp_path/'source';src.mkdir();(src/'_invocation').write_bytes(b'call frozen');(src/'_mrosource').write_bytes(b'frozen definitions')
    stage=src/'stage';stage.mkdir();chunk=stage/'chnk0-u0123456789';chunk.mkdir();(chunk/'files').mkdir();(chunk/'files/payload').write_bytes(b'not copied')
    (chunk/'_complete').write_bytes(b'complete');(chunk/'_outs').write_bytes(b'{"path":"unchanged/original/payload"}')
    try:(stage/'chnk0').symlink_to('chnk0-u0123456789',target_is_directory=True)
    except OSError:pytest.skip('Native symlinks required')
    return src,chunk

def test_alias_flattening_preserves_metadata_not_payload(tmp_path):
    src,chunk=source(tmp_path);dst=tmp_path/'copy';r=module.clone(src,dst)
    assert (dst/'stage/chnk0/_outs').read_bytes()==(chunk/'_outs').read_bytes()
    assert not (dst/'stage/chnk0').is_symlink()
    assert not (dst/'stage/chnk0-u0123456789').exists()
    assert not list((dst/'stage/chnk0/files').iterdir())
    assert (chunk/'files/payload').read_bytes()==b'not copied'
    assert r['execution_allowed'] is False and r['output_paths_rewritten'] is False

def test_active_source_refused(tmp_path):
    src,_=source(tmp_path);(src/'_lock').write_bytes(b'active')
    with pytest.raises(ValueError):module.clone(src,tmp_path/'copy')

def test_source_nested_target_refused(tmp_path):
    src,_=source(tmp_path)
    with pytest.raises(ValueError):module.clone(src,src/'copy')

def test_escaping_metadata_refused(tmp_path):
    src,_=source(tmp_path);outside=tmp_path/'outside';outside.write_bytes(b'outside')
    (src/'_external').symlink_to(outside)
    with pytest.raises(ValueError,match='Escaping metadata'):module.clone(src,tmp_path/'copy')

def test_fifo_metadata_refused(tmp_path):
    import os
    if not hasattr(os,'mkfifo'):pytest.skip('POSIX fifo only')
    src,_=source(tmp_path);os.mkfifo(src/'_pipe')
    with pytest.raises(ValueError,match='Non-regular'):module.clone(src,tmp_path/'copy')
