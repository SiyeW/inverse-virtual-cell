import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from migrate_verified_tree import plan,hash_file,inventory,require_source_scope,require_verified_target,cleanup_group


def fixture_group(tmp_path):
    source_base=tmp_path/'source'
    source=source_base/'delivery'
    source.mkdir(parents=True)
    (source/'a.bin').write_bytes(b'abc')
    target=tmp_path/'target'
    target.mkdir()
    config={'source_base':str(source_base),'target_base':str(target),
            'groups':[{'name':'test','source':str(source),'relative_target':'lab/delivery'}]}
    return config,plan(config)[0]


def test_hash_and_scope_rejects_extra_and_changes(tmp_path):
    _,g=fixture_group(tmp_path)
    source=Path(g['source'])
    digest=hash_file(source/'a.bin')
    assert digest['sha256']=='ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'
    assert digest['md5']=='900150983cd24fb0d6963f7d28e17f72'
    require_source_scope(g)
    (source/'extra').write_bytes(b'x')
    with pytest.raises(RuntimeError): require_source_scope(g)


def test_no_cleanup_without_destination_proof(tmp_path):
    config,g=fixture_group(tmp_path)
    Path(g['target']).mkdir(parents=True)
    with pytest.raises(RuntimeError):
        cleanup_group(g,Path(config['source_base']),lambda:None)
    assert (Path(g['source'])/'a.bin').read_bytes()==b'abc'


def test_destination_exact_hash_proof_and_changed_file_rejected(tmp_path):
    _,g=fixture_group(tmp_path)
    src,dst=Path(g['source']),Path(g['target'])
    dst.mkdir(parents=True)
    (dst/'a.bin').write_bytes(b'abc')
    g['source_hashes']['a.bin']=hash_file(src/'a.bin')
    g['verified']['a.bin']=hash_file(dst/'a.bin')
    require_verified_target(g)
    (dst/'a.bin').write_bytes(b'xyz')
    with pytest.raises(RuntimeError): require_verified_target(g)


def test_parent_target_and_source_overlap_rejected(tmp_path):
    config,_=fixture_group(tmp_path)
    config['groups'][0]['relative_target']='../escape'
    with pytest.raises(ValueError): plan(config)
    config['groups'][0]['relative_target']='safe'
    config['groups'].append(dict(config['groups'][0],name='duplicate'))
    with pytest.raises(ValueError): plan(config)
