"""Restrict this process tree to the standard read/write copy fallback.

On Linux x86_64, deny only copy_file_range with ENOSYS. Rust/Go then use
their existing non-offloaded copy implementation. No mount or privilege change.
"""
import ctypes, errno, os, platform, sys

class Filter(ctypes.Structure):
    _fields_=[('code',ctypes.c_ushort),('jt',ctypes.c_ubyte),('jf',ctypes.c_ubyte),('k',ctypes.c_uint32)]

class Program(ctypes.Structure):
    _fields_=[('length',ctypes.c_ushort),('filters',ctypes.POINTER(Filter))]

def enable():
    if sys.platform!='linux' or platform.machine()!='x86_64':
        raise RuntimeError('This reviewed filter is Linux x86_64 only')
    # seccomp_data: nr at offset 0, arch at offset 4. Non-x86_64 calls
    # are untouched; the x86_64 copy_file_range syscall number is 326.
    entries=(Filter*6)(
        Filter(0x20,0,0,4),                 # LD W ABS arch
        Filter(0x15,0,3,0xc000003e),       # JEQ AUDIT_ARCH_X86_64
        Filter(0x20,0,0,0),                 # LD W ABS nr
        Filter(0x15,0,1,326),              # JEQ copy_file_range
        Filter(0x06,0,0,0x00050000|errno.ENOSYS),
        Filter(0x06,0,0,0x7fff0000))        # ALLOW everything else
    program=Program(len(entries),entries)
    libc=ctypes.CDLL(None,use_errno=True)
    if libc.prctl(38,1,0,0,0)!=0:
        raise OSError(ctypes.get_errno(),'PR_SET_NO_NEW_PRIVS failed')
    if libc.prctl(22,2,ctypes.byref(program),0,0)!=0:
        raise OSError(ctypes.get_errno(),'PR_SET_SECCOMP failed')

def main():
    args=sys.argv[1:]
    if args and args[0]=='--':args=args[1:]
    if not args:raise SystemExit('Supply the owned command to execute')
    enable()
    os.execvpe(args[0],args,os.environ)

if __name__=='__main__':main()
