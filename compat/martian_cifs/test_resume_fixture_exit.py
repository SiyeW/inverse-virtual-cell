import os
import subprocess
import sys
import pytest
from run_resume_fixture import still_live,owned_tree


def test_only_original_live_writers_block_resume():
    known=[dict(pid=x,uid=2601,startticks='123',state='S') for x in range(1,5)]
    current={
        1:dict(pid=1,uid=2601,startticks='123',state='Z'),
        2:dict(pid=2,uid=2601,startticks='999',state='R'),
        3:dict(pid=3,uid=2301,startticks='123',state='R'),
        4:dict(pid=4,uid=2601,startticks='123',state='D'),
    }
    assert [x['pid'] for x in still_live(known,current.get)]==[4]
    current.pop(4)
    assert still_live(known,current.get)==[]


@pytest.mark.skipif(sys.platform!='linux',reason='Linux /proc ownership integration')
def test_known_real_process_exit_is_verified():
    proc=subprocess.Popen([sys.executable,'-c','import time;time.sleep(10)'])
    try:
        known=owned_tree(proc.pid)
        assert any(x['pid']==proc.pid and x['uid']==os.getuid() for x in known)
        assert any(x['pid']==proc.pid for x in still_live(known))
    finally:
        proc.terminate();proc.wait(timeout=10)
    assert not still_live(known)
