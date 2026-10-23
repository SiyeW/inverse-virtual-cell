"""Replay accepted, bounded RNA+guide prefixes; never regenerate/read raw FASTQs."""
import argparse,csv,hashlib,json,os,shutil,signal,subprocess,sys,time
from pathlib import Path
from verify_fixture import compare_matrix,inventory
import copy_fallback

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def validate_recipe(status,receipt,rows):
    if status.get('stage')!='real_library_prefix_engineering_passed' or status.get('scientific_result') is not False:
        raise ValueError('Accepted engineering prefix required')
    if status.get('forced_cells')!=10 or receipt.get('pairs_per_library')!=100000 or receipt.get('scientific_result') is not False:
        raise ValueError('Only the frozen 100k-pair/10-cell engineering recipe is admitted')
    libs=receipt['libraries']
    if len(libs)!=2 or {x['feature_type'] for x in libs}!={'Gene Expression','CRISPR Guide Capture'}:
        raise ValueError('Expected paired RNA and CRISPR libraries')
    expected=[dict(fastqs=x['fastq_dir'],sample=x['fastq_id'],library_type=x['feature_type']) for x in libs]
    if rows!=expected:raise ValueError('Libraries CSV does not point exclusively to accepted prefixes')

def run(baseline,target,tool):
    baseline=baseline.resolve(strict=True);target=target.absolute();tool=tool.resolve(strict=True)
    if target.exists() or target.is_symlink() or 'provenance' not in target.parts:
        raise ValueError('Unused project provenance target required')
    view=json.loads((tool/'compat_view.json').read_text())
    for item in view['files']:
        if item['kind'] in ['verified_copy','fixed_loader_config'] and digest(tool/item['path'])!=item['sha256']:
            raise ValueError('Frozen compatibility software differs')
    status=json.loads((baseline/'chip02_prefix_status.json').read_text())
    receipt=json.loads((baseline/'chip02_prefix_100k/prefix_receipt.json').read_text())
    argv=list(status['argv'])
    with Path(argv[argv.index('--libraries')+1]).open(newline='') as f:rows=list(csv.DictReader(f))
    validate_recipe(status,receipt,rows)
    if argv[argv.index('--force-cells')+1]!='10' or '--create-bam=false' not in argv or '--localcores=2' not in argv or '--localmem=48' not in argv:
        raise ValueError('Prefix execution recipe differs')
    total=0
    for lib in receipt['libraries']:
        if len(lib['files'])!=2:raise ValueError('Expected paired prefix files')
        for item in lib['files']:
            p=Path(item['output']).resolve(strict=True)
            if not p.is_relative_to(baseline/'chip02_prefix_100k'):raise ValueError('Not an owned prefix')
            total+=p.stat().st_size
            if total>50_000_000 or p.stat().st_size!=item['output_bytes'] or digest(p)!=item['output_sha256']:
                raise ValueError('Prefix hash/size/budget mismatch')
    ref=Path(argv[argv.index('--transcriptome')+1]);proof=json.loads((baseline/'reference_content_receipt.json').read_text())
    if proof.get('complete') is not True or str(ref)!=proof['reference_root']:raise ValueError('Reference proof missing')
    for rel,item in proof['files'].items():
        s=(ref/rel).stat()
        if (s.st_size,s.st_mtime_ns)!=(item['bytes'],item['mtime_ns']):raise ValueError('Reference identity changed')
    accepted=json.loads((baseline/'chip02_prefix_acceptance.json').read_text())
    old=baseline/'chip02_count_100k/outs'
    old_files=inventory(old)
    if {x['path']:{k:v for k,v in x.items() if k!='path'} for x in old_files}!=accepted['outputs']:
        raise ValueError('Accepted baseline outputs changed')
    feature=Path(argv[argv.index('--feature-ref')+1])
    if digest(feature)!=accepted['outputs']['feature_reference.csv']['sha256']:raise ValueError('Guide reference changed')
    available=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024
    if available<64*1024**3:raise ValueError('Insufficient shared memory headroom')
    target.mkdir();output=target/'count';argv[0]=str(tool/'cellranger');argv[argv.index('--output-dir')+1]=str(output)
    env=dict(os.environ,MRO_UNIQUIFIED_DIRECTORIES='disable',MRO_CIFS_COMPAT='1',MRO_CIFS_SAME_MODE_ROOT=str(target))
    command=[sys.executable,str(Path(copy_fallback.__file__)),'--']+argv
    rec=dict(stage='running_prefix_count',pid=os.getpid(),uid=os.getuid(),argv=argv,scientific_result=False,original_fastq_bytes_read=0,prefix_bytes_verified=total,forced_cells=10,
             helper_sha256=digest(Path(__file__)),software_view_sha256=digest(tool/'compat_view.json'),
             prefix_receipt_sha256=digest(baseline/'chip02_prefix_100k/prefix_receipt.json'),
             reference_receipt_sha256=digest(baseline/'reference_content_receipt.json'))
    start=time.monotonic();proc=None
    with (target/'owner.lock').open('x') as f:json.dump(rec,f)
    def save():
        rec['elapsed_seconds']=time.monotonic()-start
        temp=target/'status.tmp';temp.write_text(json.dumps(rec,indent=2)+'\n');temp.replace(target/'status.json')
    try:
        with (target/'count.stdout').open('w') as out,(target/'count.stderr').open('w') as err:
            proc=subprocess.Popen(command,cwd=target,env=env,stdout=out,stderr=err,start_new_session=True)
            rec['pipeline_pid']=proc.pid;save()
            while proc.poll() is None:
                if time.monotonic()-start>1800 or shutil.disk_usage(target).free<20_000_000_000:raise RuntimeError('Test budget reached')
                time.sleep(5);save()
        if proc.returncode:raise RuntimeError('Prefix replay failed; preserve without automatic retry')
        if 'Could not move output files' in (output/'_log').read_text():raise ValueError('Publication error')
        json.loads((output/'SC_RNA_COUNTER_CS/fork0/_outs').read_text())
        new=output/'outs';new_files=inventory(new)
        if [x['path'] for x in old_files]!=[x['path'] for x in new_files]:raise ValueError('Output coverage differs')
        comparisons={n:compare_matrix(old/n,new/n) for n in ['raw_feature_bc_matrix.h5','filtered_feature_bc_matrix.h5']}
        for n in ['crispr_analysis/protospacer_calls_per_cell.csv','feature_reference.csv','metrics_summary.csv']:
            if (old/n).read_bytes()!=(new/n).read_bytes():raise ValueError('Guide/QC identity differs: '+n)
        result=dict(complete=True,scientific_result=False,original_fastq_bytes_read=0,matrices=comparisons,guide_calls_equal=True,feature_reference_equal=True,metrics_equal=True,outputs=new_files,no_bam=True)
        (target/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
        rec['stage']='prefix_compatibility_passed'
    except Exception as exc:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid,signal.SIGINT)
            try:proc.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=15)
        rec.update(stage='failed_preserved',error=repr(exc));raise
    finally:
        save();(target/'owner.lock').unlink()
    print(json.dumps(rec))

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--target',type=Path,required=True);p.add_argument('--tool',type=Path,required=True)
    a=p.parse_args();run(a.baseline,a.target,a.tool)

if __name__=='__main__':main()
