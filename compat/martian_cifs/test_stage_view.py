import importlib.util
from pathlib import Path
import pytest
import json
import hashlib

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

def test_report_overlay_is_independent(tmp_path):
    src,runtime,shim=fixture(tmp_path);report=src/module.REPORT_PATH;report.parent.mkdir(parents=True)
    source=b'import shutil\n\ndef main(args, outs):\n    shutil.copy(args.feature_reference, outs.feature_reference)\n'
    report.write_bytes(source);helper=tmp_path/'helper.py';helper.write_text('def checked_copy(a,b): pass\n')
    out=tmp_path/'view';module.stage(src,runtime,out,shim,helper,hashlib.sha256(source).hexdigest())
    assert report.read_bytes()==source
    assert b'cifs_content_copy.checked_copy' in (out/module.REPORT_PATH).read_bytes()
    assert not (out/module.REPORT_PATH).is_symlink()
    assert json.loads((out/'.env.json').read_text())['envs'][-1]=={'key':'PYTHONPATH','type':'path_prepend','value':'compat'}

def test_directory_overlay_preserves_source_and_scientific_callback(tmp_path):
    src,runtime,shim=fixture(tmp_path)
    target=src/module.CR_IO_PATH;target.parent.mkdir(parents=True)
    source=b'import shutil\n\ndef hardlink_with_fallback(src, dst):\n    def _hardlink_file_with_fallback(a,b): pass\n    shutil.copytree(src, dst, copy_function=_hardlink_file_with_fallback, dirs_exist_ok=True)\n'
    target.write_bytes(source);helper=tmp_path/'directory.py';helper.write_text('def copytree(a,b,**k): pass\n')
    out=tmp_path/'view'
    record=module.stage(src,runtime,out,shim,directory_helper=helper,cr_io_sha256=hashlib.sha256(source).hexdigest())
    assert target.read_bytes()==source
    assert b'cifs_directory_copy.copytree' in (out/module.CR_IO_PATH).read_bytes()
    assert not (out/module.CR_IO_PATH).is_symlink()
    assert record['scientific_binaries_modified'] is False
    assert record['directory_publication_policy']=='content_verified_inherit_destination'
