import ctypes,errno,os,platform,subprocess,sys
from pathlib import Path
import pytest

@pytest.mark.skipif(sys.platform!='linux' or platform.machine()!='x86_64',reason='Reviewed Linux x86_64 filter only')
def test_filter_changes_child_only_and_preserves_regular_io(tmp_path):
    libc=ctypes.CDLL(None,use_errno=True)
    libc.syscall(326,-1,0,-1,0,1,0);before=ctypes.get_errno()
    code='''import ctypes,errno,pathlib,sys
c=ctypes.CDLL(None,use_errno=True)
assert c.syscall(326,-1,0,-1,0,1,0)==-1
assert ctypes.get_errno()==errno.ENOSYS
p=pathlib.Path(sys.argv[1]);p.write_bytes(b"ordinary IO preserved")
assert p.read_bytes()==b"ordinary IO preserved"
print("child fallback active; ordinary IO preserved")
'''
    subprocess.run([sys.executable,str(Path(__file__).with_name('copy_fallback.py')),'--',sys.executable,'-c',code,str(tmp_path/'payload')],check=True,timeout=10)
    libc.syscall(326,-1,0,-1,0,1,0)
    assert ctypes.get_errno()==before,'Parent process was changed'
