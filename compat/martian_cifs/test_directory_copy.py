import hashlib
import ast
import errno
import json
import os
from pathlib import Path
import shutil
import typing
import pytest
import cifs_directory_copy as copy
from adapt_directory_copy import adapt


def setup_tree(tmp_path, monkeypatch):
    monkeypatch.setenv('MRO_CIFS_COMPAT', '1')
    monkeypatch.setenv('MRO_CIFS_SAME_MODE_ROOT', str(tmp_path))
    monkeypatch.setattr(copy, '_cifs', lambda root: True)
    src = tmp_path / 'source';src.mkdir()
    (src/'nested').mkdir();(src/'nested'/'counts.csv').write_bytes(b'RNA,UMI\nA,42\n')
    return src, tmp_path/'destination'


def test_verified_content_and_no_directory_metadata(tmp_path, monkeypatch):
    src, dst = setup_tree(tmp_path, monkeypatch)
    monkeypatch.setattr(shutil, 'copystat', lambda *a, **k: (_ for _ in ()).throw(PermissionError(1, 'metadata denied')))
    copied=[]
    def callback(a,b):
        copied.append(a);return shutil.copyfile(a,b)
    copy.copytree(bytes(src), bytes(dst), copy_function=callback, dirs_exist_ok=True)
    assert (dst/'nested/counts.csv').read_bytes()==b'RNA,UMI\nA,42\n'
    assert (src/'nested/counts.csv').read_bytes()==b'RNA,UMI\nA,42\n'
    copy.copytree(src,dst,copy_function=callback,dirs_exist_ok=True)
    assert len(copied)==1


def test_bad_existing_content_and_failed_copy_are_rejected(tmp_path, monkeypatch):
    src,dst=setup_tree(tmp_path,monkeypatch)
    dst.mkdir();(dst/'nested').mkdir();(dst/'nested/counts.csv').write_bytes(b'wrong')
    with pytest.raises(ValueError,match='Existing publication'):
        copy.copytree(src,dst,copy_function=shutil.copyfile,dirs_exist_ok=True)
    assert (dst/'nested/counts.csv').read_bytes()==b'wrong'
    with pytest.raises(PermissionError):
        copy.copytree(src,tmp_path/'denied',copy_function=lambda *a: (_ for _ in ()).throw(PermissionError(13,'data denied')),dirs_exist_ok=True)


def test_escape_and_links_are_rejected(tmp_path,monkeypatch):
    src,dst=setup_tree(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='escapes'):
        copy.copytree(src,tmp_path.parent/'outside',copy_function=shutil.copyfile,dirs_exist_ok=True)
    try:(src/'link').symlink_to(src/'nested/counts.csv')
    except OSError:pytest.skip('Symlink privilege unavailable')
    with pytest.raises(ValueError,match='escapes'):
        copy.copytree(src,dst,copy_function=shutil.copyfile,dirs_exist_ok=True)


def test_no_opt_in_or_native_filesystem_keeps_standard_behavior(tmp_path,monkeypatch):
    src,dst=setup_tree(tmp_path,monkeypatch)
    calls=[]
    monkeypatch.setattr(shutil,'copytree',lambda *a,**k:calls.append((a,k)))
    monkeypatch.delenv('MRO_CIFS_COMPAT')
    copy.copytree(src,dst,copy_function=shutil.copyfile,dirs_exist_ok=True)
    monkeypatch.setenv('MRO_CIFS_COMPAT','1');monkeypatch.setattr(copy,'_cifs',lambda root:False)
    copy.copytree(src,dst,copy_function=shutil.copyfile,dirs_exist_ok=True)
    assert len(calls)==2


def test_pinned_single_call_adaptation():
    source=b'import shutil\n\ndef hardlink_with_fallback(src, dst):\n    def _hardlink_file_with_fallback(a,b): pass\n    shutil.copytree(src, dst, copy_function=_hardlink_file_with_fallback, dirs_exist_ok=True)\n'
    h=hashlib.sha256(source).hexdigest()
    result=adapt(source,h)
    assert b'cifs_directory_copy.copytree' in result
    with pytest.raises(ValueError):adapt(source,'0'*64)
    with pytest.raises(ValueError):adapt(result,hashlib.sha256(result).hexdigest())


@pytest.mark.skipif(not os.environ.get('MRO_CIFS_DIRECTORY_TEST_ROOT'),reason='Explicit real CIFS fixture only')
def test_actual_vendor_directory_publication(monkeypatch):
    root=Path(os.environ['MRO_CIFS_DIRECTORY_TEST_ROOT'])
    root.mkdir(exist_ok=False)
    monkeypatch.setenv('MRO_CIFS_COMPAT','1')
    monkeypatch.setenv('MRO_CIFS_SAME_MODE_ROOT',str(root))
    assert copy._cifs(root)
    source=Path(os.environ['MRO_CIFS_CR_IO_SOURCE']).read_bytes()
    before=hashlib.sha256(source).hexdigest()
    expected=os.environ['MRO_CIFS_CR_IO_SHA256'];assert before==expected
    def load(raw):
        node=next(n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name=='hardlink_with_fallback')
        space=dict(os=os,shutil=shutil,errno=errno,AnyStr=typing.AnyStr,cifs_directory_copy=copy)
        exec(compile(ast.Module(body=[node],type_ignores=[]),'<pinned vendor hardlink>','exec'),space)
        return space['hardlink_with_fallback']
    src=root/'source';src.mkdir();(src/'nested').mkdir()
    (src/'nested'/'counts.csv').write_bytes(b'gene,UMI\nRNA1,42\n')
    original_error=None
    try:load(source)(str(src),str(root/'original_destination'))
    except shutil.Error as exc:original_error=str(exc)
    # Record rather than invent the host's current metadata behavior.
    corrected=load(adapt(source,expected))
    corrected(os.fsencode(src),os.fsencode(root/'fixed_destination'))
    assert (root/'fixed_destination/nested/counts.csv').read_bytes()==(src/'nested/counts.csv').read_bytes()
    corrected(str(src),str(root/'fixed_destination'))
    (root/'receipt.json').write_text(json.dumps(dict(original_error=original_error,content_equal=True,repeat_equal=True,source_sha256=before)))
