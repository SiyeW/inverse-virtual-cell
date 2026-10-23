"""Run only against an unused, owned test root on the actual CIFS mount."""
import argparse, hashlib, json, os, pathlib, subprocess, sys

def child(path, mode, expected):
    p=pathlib.Path(path);before=p.read_bytes();old=p.stat().st_mode&0o7777
    with p.open('rb') as f:
        try: os.fchmod(f.fileno(),int(mode,8));actual='success'
        except PermissionError: actual='denied'
    assert actual==expected,(actual,expected)
    assert p.read_bytes()==before and p.stat().st_mode&0o7777==old
    print(json.dumps(dict(path=str(p),mode=mode,result=actual,content_sha256=hashlib.sha256(before).hexdigest())))

def main():
    if len(sys.argv)>1 and sys.argv[1]=='child':
        child(*sys.argv[2:]);return
    a=argparse.ArgumentParser();a.add_argument('--root',type=pathlib.Path,required=True);a.add_argument('--library',type=pathlib.Path,required=True);args=a.parse_args()
    args.root.mkdir();inside=args.root/'inside';outside=args.root/'outside';inside.mkdir();outside.mkdir()
    p=inside/'payload';q=outside/'payload';p.write_bytes(b'unchanged-test-content');q.write_bytes(b'outside-content')
    assert p.stat().st_mode&0o7777==0o777,'Test requires observed CIFS mode 0777'
    clean=dict(os.environ);clean.pop('LD_PRELOAD',None);clean.pop('MRO_CIFS_SAME_MODE_ROOT',None)
    cases=[(p,'777','denied',clean),(p,'777','success',dict(clean,LD_PRELOAD=str(args.library),MRO_CIFS_SAME_MODE_ROOT=str(inside))),
           (p,'700','denied',dict(clean,LD_PRELOAD=str(args.library),MRO_CIFS_SAME_MODE_ROOT=str(inside))),
           (q,'777','denied',dict(clean,LD_PRELOAD=str(args.library),MRO_CIFS_SAME_MODE_ROOT=str(inside))),
           (p,'777','denied',dict(clean,LD_PRELOAD=str(args.library)))]
    for path,mode,expected,env in cases:
        subprocess.run([sys.executable,__file__,'child',str(path),mode,expected],env=env,check=True,timeout=15)
    print('5 scoped fchmod tests passed; file bytes and modes unchanged')

if __name__=='__main__':main()
