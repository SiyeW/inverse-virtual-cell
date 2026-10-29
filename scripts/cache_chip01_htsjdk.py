"""Read an existing local indexed BAM, then transfer only a small RNA cache.

Designed to run in an independent Windows host; no complete BAM upload, new
BAM, mount/network change or inference of DNA/functional editing efficiency.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,socket,struct,subprocess,sys,time,zipfile

from cache_guided_rna_alignments import canonical_hash,sha_file,atomic,frozen_identity,inventory


def sid():
    return subprocess.check_output(['powershell','-NoProfile','-Command',
        '[System.Security.Principal.WindowsIdentity]::GetCurrent().User.Value'],text=True).strip()


def source_header(path):
    with gzip.open(path,'rb') as f:
        if f.read(4)!=b'BAM\x01':raise ValueError('Source is not BAM')
        n=struct.unpack('<i',f.read(4))[0]
        if not 0<n<131072:raise ValueError('Unbounded source header')
        text=f.read(n).decode();count=struct.unpack('<i',f.read(4))[0];refs=[]
        for _ in range(count):
            size=struct.unpack('<i',f.read(4))[0]
            name=f.read(size).rstrip(b'\0').decode();length=struct.unpack('<i',f.read(4))[0]
            refs.append(dict(name=name,length=length))
    return text,refs


def verify_source(task):
    if sid()!=task['source_owner_sid'] or socket.gethostname()!=task['source_host']:raise ValueError('Windows source owner/host differs')
    for name,expected in [('source_BAM','source_BAM_identity'),('source_index','source_index_identity')]:
        if inventory(task[name])!=task[expected]:raise ValueError('Local source identity differs')
    text,refs=source_header(task['source_BAM'])
    if refs!=task['reference_contigs_ordered'] or hashlib.sha256(text.encode()).hexdigest()!=task['source_BAM_raw_header_sha256']:
        raise ValueError('Complete ordered local reference identity differs')
    return text


def verified_ready(path,staging,by_region,identity):
    item=json.loads(path.read_text());region=item['region'];data=staging/(region['region_id']+'.jsonl.gz')
    if region!=by_region[region['region_id']]:raise ValueError('Old region is not equivalent')
    if data.stat().st_size!=item['bytes'] or sha_file(data)!=item['sha256']:raise ValueError('Ready cache content differs')
    with gzip.open(data,'rt') as f:header=json.loads(next(f))
    if header!=dict(kind='RNA_SAM_cache_header',identity=identity,region=region):raise ValueError('Old cache source differs')
    return item


def main():
    p=argparse.ArgumentParser();p.add_argument('--job',type=Path,required=True);args=p.parse_args()
    cfg=json.loads(args.job.read_text());task_path=Path(cfg['task_manifest']);task=json.loads(task_path.read_text());op=Path(cfg['operation_root'])
    lock=op/'chip01_source.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    import psutil
    os.write(fd,json.dumps(dict(task_id=task['task_id'],dataset_id='chip01',pid=os.getpid(),
        creation_time=psutil.Process().create_time(),owner_sid=sid(),argv=sys.argv)).encode());os.close(fd)
    def state(stage,**kw):atomic(op/'chip01_job_status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kw))
    try:
        for path,sha in cfg['immutable_files'].items():
            if sha_file(path)!=sha:raise ValueError('Fixed source/metadata file differs: '+path)
        text=verify_source(task)
        candidate=json.loads(Path(cfg['candidate_contract']).read_text());identity=frozen_identity(task,candidate)
        staging=Path(cfg['staging']);staging.mkdir(exist_ok=True)
        (staging/'source_header.sam').write_bytes(text.encode('utf8'))
        receipts=[];ready_files=list(staging.glob('*.ready.json'));skip=[]
        by_region={r['region_id']:r for r in candidate['regions']}
        for path in ready_files:
            item=verified_ready(path,staging,by_region,identity)
            skip.append(item['region']['region_id'])
        java_cfg=dict(mode='extract',bam=task['source_BAM'],bam_index=task['source_index'],
            bam_identity=task['source_BAM_identity'],bam_index_identity=task['source_index_identity'],
            reference_contigs_ordered=task['reference_contigs_ordered'],regions=candidate['regions'],
            identity=identity,assignments_tsv=cfg['assignments_tsv'],NTC_guides_txt=cfg['NTC_guides_txt'],
            local_staging=str(staging),max_fetched_records_per_region=200000,max_local_staging_bytes=32*1024**2,
            wait_for_upload=False,skip_regions=skip)
        atomic(op/'chip01_java_config.json',java_cfg)
        command=[cfg['java'],'-Xmx256m','-XX:ActiveProcessorCount=2','--class-path',cfg['classpath'],cfg['java_source'],str(op/'chip01_java_config.json')]
        state('local_indexed_RNA_extraction',completed_regions=len(skip),total_regions=len(by_region))
        with (op/'chip01_java.stdout').open('w') as out,(op/'chip01_java.stderr').open('w') as err:
            result=subprocess.run(command,stdout=out,stderr=err,timeout=3600,creationflags=0x08000000)
        if result.returncode:raise ValueError('HTSJDK extraction failed; preserve all ready caches')
        verify_source(task)
        files={};total=0
        for region in candidate['regions']:
            rid=region['region_id'];ready=verified_ready(staging/(rid+'.ready.json'),staging,by_region,identity);path=staging/(rid+'.jsonl.gz')
            name='cache/'+path.name;proof=dict(bytes=path.stat().st_size,sha256=sha_file(path));files[name]=proof;total+=proof['bytes']
            receipts.append(dict(identity_sha256=canonical_hash(identity),region=region,file=name,**proof,
                                 fetched_records=ready['fetched_records'],selected_records=ready['selected_records'],
                                 selected_without_UB=ready['selected_without_UB'],scan_capped=ready['scan_capped'],seconds=ready['seconds']))
        if total>32*1024**2:raise ValueError('Small local staging bound exceeded; retained for scoped transfer recovery')
        files['source_header.sam']=dict(bytes=(staging/'source_header.sam').stat().st_size,sha256=sha_file(staging/'source_header.sam'))
        archive=op/'chip01_cache_import.zip'
        with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_STORED) as z:
            for name in files:z.write(staging/(Path(name).name),name)
        manifest=dict(task_id=task['task_id'],dataset_id='chip01',source_task_manifest_sha256=sha_file(task_path),
            identity=identity,source_host=task['source_host'],source_owner_sid=task['source_owner_sid'],
            source_final_identity_verified=True,archive_sha256=sha_file(archive),files=files,regions=receipts,
            raw_BAM_copied=False,new_BAM_created=False)
        atomic(op/'chip01_import_manifest.json',manifest)
        state('upload_small_verified_cache',regions=len(receipts),cache_bytes=total)
        native=task['native_directory'];names=[archive,op/'chip01_import_manifest.json']
        batch=''.join('put "'+p.as_posix()+'" "'+native+'/'+p.name+'"\n' for p in names)
        c=subprocess.run(['sftp','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','ConnectionAttempts=1','-b','-','rhesusbase'],input=batch,text=True,capture_output=True,timeout=180,creationflags=0x08000000)
        if c.returncode:raise ValueError('Cache upload failed; source/cache retained: '+c.stderr[-500:])
        import shlex
        command=[cfg['server_python'],cfg['server_receiver'],'--task',native+'/task_manifest.json','--candidates',native+'/candidate_contract.json','--manifest',native+'/chip01_import_manifest.json','--archive',native+'/chip01_cache_import.zip']
        c=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','ConnectionAttempts=1','rhesusbase',shlex.join(command)],capture_output=True,text=True,timeout=180,creationflags=0x08000000)
        if c.returncode:raise ValueError('Import rejected; raw and derived cache remain separate: '+c.stderr[-1200:])
        (op/'chip01_import.stdout').write_text(c.stdout);state('source_cache_import_complete',regions=len(receipts),cache_bytes=total,server_result=json.loads(c.stdout))
    except BaseException as e:state('failed_pending_recovery',error=repr(e));raise
    finally:lock.unlink()


if __name__=='__main__':main()
