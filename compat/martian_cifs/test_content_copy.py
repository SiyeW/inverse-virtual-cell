import hashlib,importlib.util,sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).parent))
import cifs_content_copy,adapt_report

def test_copy_scope_hash_and_source_preservation(tmp_path,monkeypatch):
    root=tmp_path/'out';root.mkdir();src=tmp_path/'source';src.write_bytes(b'guide identity retained')
    old=src.stat();monkeypatch.setenv('MRO_CIFS_COMPAT','1');monkeypatch.setenv('MRO_CIFS_SAME_MODE_ROOT',str(root))
    dst=root/'feature.csv';cifs_content_copy.checked_copy(src,dst);cifs_content_copy.checked_copy(src,dst)
    assert src.read_bytes()==dst.read_bytes()==b'guide identity retained'
    assert src.stat().st_mode==old.st_mode
    with pytest.raises(ValueError):cifs_content_copy.checked_copy(src,tmp_path/'outside')
    dst.write_bytes(b'collision')
    with pytest.raises(ValueError):cifs_content_copy.checked_copy(src,dst)

def test_report_transform_is_narrow():
    source=b'import shutil\n\ndef main(args, outs):\n    outs.count = args.count + 1\n    shutil.copy(args.feature_reference, outs.feature_reference)\n'
    h=hashlib.sha256(source).hexdigest();result=adapt_report.adapt(source,h)
    assert b'outs.count = args.count + 1' in result
    assert b'cifs_content_copy.checked_copy' in result
    with pytest.raises(ValueError):adapt_report.adapt(source,'0'*64)
    with pytest.raises(ValueError):adapt_report.adapt(source.replace(b'shutil.copy',b'shutil.copyfile'),hashlib.sha256(source.replace(b'shutil.copy',b'shutil.copyfile')).hexdigest())
