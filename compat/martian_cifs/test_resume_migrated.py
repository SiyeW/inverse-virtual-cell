from pathlib import Path
import pytest
from resume_migrated import check_inputs,unexpected_completed_reexecution
from migrate_checkpoint import identity

def test_identity_allows_only_inode_drift(tmp_path):
    p=tmp_path/'input';p.write_text('original');s=identity(p);s['inode']+=100
    check_inputs([dict(path=str(p),identity=s)])
    p.write_text('changed content')
    with pytest.raises(ValueError):check_inputs([dict(path=str(p),identity=s)])

def test_reexecution_guard_distinguishes_chunk_and_history():
    prefix='[runtime] '
    assert unexpected_completed_reexecution(prefix+'(run:local) ID.chip.SC_MULTI.MAKE_SHARD.fork0.chnk1.main')
    assert not unexpected_completed_reexecution(prefix+'(run:local) ID.chip.SC_MULTI.MAKE_SHARD.fork0.chnk0.main')
    assert not unexpected_completed_reexecution(prefix+'(reset-partial) ID.chip.SC_MULTI.MAKE_SHARD.fork0.chnk1')
