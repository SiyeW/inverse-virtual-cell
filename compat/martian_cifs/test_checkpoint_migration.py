import json
from pathlib import Path
import pytest
import migrate_checkpoint as m

def setup(tmp_path):
    src=tmp_path/'old';src.mkdir();(src/'_invocation').write_text('frozen');(src/'_mrosource').write_text('frozen')
    chunk=src/'stage/chnk1-u0123456789';(chunk/'files').mkdir(parents=True)
    try:(src/'stage/chnk1').symlink_to(chunk.name,target_is_directory=True)
    except OSError:pytest.skip('Native symlink required')
    payload=chunk/'files/done';payload.write_bytes(b'done-result')
    (chunk/'_complete').write_text('done');(chunk/'_outs').write_text(json.dumps({'result':str(payload).replace('\\','/')}))
    (src/'failed/files').mkdir(parents=True);(src/'failed/files/partial').write_bytes(b'do not copy')
    plan=dict(source=src.resolve().as_posix(),target=(tmp_path/'new').as_posix(),expected_files=1,expected_bytes=payload.stat().st_size,
              files=[dict(path=payload.as_posix(),identity=m.identity(payload))])
    return src,payload,plan

@pytest.mark.skipif(__import__('os').name!='posix',reason='Linux execution paths')
def test_complete_only_and_metadata_rebased(tmp_path):
    src,p,plan=setup(tmp_path);before=m.sha(p);r=m.migrate(plan,tmp_path/'control')
    dst=Path(plan['target']);assert r['stage']=='migrated_pending_recovery_validation'
    assert m.sha(dst/'stage/chnk1/files/done')==before and m.sha(p)==before
    assert not list((dst/'failed/files').iterdir())
    assert json.loads((dst/'stage/chnk1/_outs').read_text())['result']==str(dst/'stage/chnk1/files/done')
    assert str(src) in (src/'stage/chnk1/_outs').read_text()

@pytest.mark.skipif(__import__('os').name!='posix',reason='Linux execution paths')
def test_changed_source_refused(tmp_path):
    src,p,plan=setup(tmp_path);p.write_bytes(b'changed')
    with pytest.raises(ValueError,match='identity'):m.migrate(plan,tmp_path/'control')
    assert not Path(plan['target']).exists()

def test_mapping_longest_first():
    maps=m.mappings(dict(source='/old',target='/new',aliases_flattened=[dict(source_target='stage/chnk1-u0123456789',path='stage/chnk1')]))
    assert m.mapped_path('/old/stage/chnk1-u0123456789/files/x',maps)=='/new/stage/chnk1/files/x'
    assert m.rebase(b'{"p":"/old/stage/chnk1-u0123456789/files/x"}',maps)==b'{"p":"/new/stage/chnk1/files/x"}'
