from pathlib import Path
import hashlib
import pytest
from resume_preserved import preserved_metadata_matches,completed_node_reexecution,owned_relative


def receipt(path,run):
    s=path.stat()
    return dict(relative=path.relative_to(run).as_posix(),bytes=s.st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),mtime_ns=s.st_mtime_ns,device=s.st_dev)


def fixture(tmp_path):
    node=tmp_path/'ALIGN_AND_COUNT/fork0/chnk00';node.mkdir(parents=True)
    marker=node/'_complete';marker.write_text('done')
    output=node/'molecules';output.write_bytes(b'protected result')
    for name in ['_uuid','_invocation']:(tmp_path/name).write_text(name)
    return dict(preserved_complete_nodes=[dict(relative=node.relative_to(tmp_path).as_posix(),complete=True,
        metadata={'_complete':receipt(marker,tmp_path)},outputs=[receipt(output,tmp_path)])],
        UUID=receipt(tmp_path/'_uuid',tmp_path),invocation=receipt(tmp_path/'_invocation',tmp_path))


def test_complete_output_and_metadata_change_are_refused(tmp_path):
    inventory=fixture(tmp_path)
    preserved_metadata_matches(tmp_path,inventory)
    (tmp_path/'ALIGN_AND_COUNT/fork0/chnk00/molecules').write_bytes(b'changed')
    with pytest.raises(ValueError,match='output identity'):preserved_metadata_matches(tmp_path,inventory)


def test_changed_completion_marker_is_refused(tmp_path):
    inventory=fixture(tmp_path)
    (tmp_path/'ALIGN_AND_COUNT/fork0/chnk00/_complete').write_text('fake')
    with pytest.raises(ValueError,match='metadata changed'):preserved_metadata_matches(tmp_path,inventory)


def test_completed_node_execution_guard_preserves_zero_padding(tmp_path):
    inventory=fixture(tmp_path)
    assert completed_node_reexecution('[runtime] (run:local) ID.chip.ALIGN_AND_COUNT.fork0.chnk00.main','chip',inventory)
    assert not completed_node_reexecution('[runtime] (run:local) ID.chip.ALIGN_AND_COUNT.fork0.chnk46.main','chip',inventory)
    assert completed_node_reexecution('[runtime] (run:local) ID.chip.MAKE_SHARD.fork0.chnk0.main','chip',inventory)
    assert completed_node_reexecution('[runtime] (run:local) ID.chip.BARCODE_CORRECTION.fork0.chnk000.main','chip',inventory)
    assert not completed_node_reexecution('[runtime] (update) ID.chip.ALIGN_AND_COUNT.fork0 (46/50 completed)','chip',inventory)


def test_preserved_path_cannot_escape_run(tmp_path):
    with pytest.raises(ValueError,match='Nonrelative'):owned_relative(tmp_path,'../outside')
    with pytest.raises(ValueError,match='Nonrelative'):owned_relative(tmp_path,str(tmp_path/'absolute'))
