"""Cache indexed RNA alignments for every frozen guide window, without new BAM.

The cache preserves original SAM records/tags. It is input to RNA allele
inspection, not an editing/KO result. Completed region receipts bind exact
source/reference/assignment inventories and the producer commit.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time


def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def atomic(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2),encoding='utf8');tmp.replace(path)


def inventory(path):
    s=Path(path).stat();return dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def frozen_identity(task,candidates):
    value=dict(task_id=task['task_id'],producer_commit=task['code_version'],
                BAM=task['source_BAM'],BAM_identity=task['source_BAM_identity'],
                BAM_index=task['source_index'],BAM_index_identity=task['source_index_identity'],
                BAM_raw_header_sha256=task['source_BAM_raw_header_sha256'],
                reference_metadata_sha256=task['reference_metadata_sha256'],
                reference_fasta_identity=task['reference_identity'],
                reference_index_identity=task['reference_index_identity'],
                reference_contigs_ordered=task['reference_contigs_ordered'],
                guide_contract_sha256=task['guide_contract_sha256'],
                assignment_sha256=task['assignment_sha256'],
                candidates_sha256=canonical_hash(candidates))
    for key in ['dataset_id','source_host','source_input_mode']:
        if key in task:value[key]=task[key]
    return value


def validate_region_receipt(item,identity,region,path):
    if item['identity_sha256']!=canonical_hash(identity) or item['region']!=region:
        raise ValueError('Cached source/reference/region/code identity differs')
    if inventory(path)['bytes']!=item['bytes'] or sha_file(path)!=item['sha256']:
        raise ValueError('Cached region content differs')
    with gzip.open(path,'rt') as f:
        header=json.loads(next(f))
    if header!=dict(kind='RNA_SAM_cache_header',identity=identity,region=region):
        raise ValueError('Cached region header inventory differs')


def check_sources(task):
    pairs=[('source_BAM','source_BAM_identity'),('source_index','source_index_identity'),
           ('reference_fasta','reference_identity'),('reference_index','reference_index_identity')]
    for path_key,identity_key in pairs:
        if inventory(task[path_key])!=task[identity_key]:
            raise ValueError(path_key+' changed; original source stays unadopted')


def samples_for_read(read,region,assigned,ntc_guides):
    if not read.has_tag('CB'):return None
    cb=read.get_tag('CB');guide=assigned.get(cb)
    if guide is None:return None
    if guide in ntc_guides:return 'NTC'
    for g in region['guides']:
        if guide==g['guide'] and read.reference_start<g['end0']+100 and read.reference_end>max(0,g['start0']-100):
            return g['sample_id']
    return None


def run(task,candidates,assignment_path,guide_path,progress):
    import pysam
    if task['task_id']!='20261030_all_target_editing':raise ValueError('Wrong parent task')
    if os.getuid()!=task['expected_owner_uid'] or socket.gethostname()!=task['verified_actual_hostname']:
        raise ValueError('Execution host or owner differs')
    if not task['code_version']:raise ValueError('No fixed pushed producer version')
    output=Path(task['canonical_output_root'])
    expected=Path('/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261030_all_target_editing_r1')
    if output!=expected:raise ValueError('Output is outside the frozen task namespace')
    check_sources(task)
    if sha_file(assignment_path)!=task['assignment_sha256']:raise ValueError('Assignment content differs')
    if sha_file(guide_path)!=task['guide_contract_sha256']:raise ValueError('Guide catalogue content differs')
    opener=gzip.open if str(assignment_path).endswith('.gz') else open
    with opener(assignment_path,'rt') as f:rows=list(csv.DictReader(f))
    assigned={r['cell_barcode']:r['guide_id'] for r in rows}
    if len(assigned)!=72107 or len(assigned)!=len(rows):raise ValueError('Wrong barcode inventory')
    with Path(guide_path).open() as f:guides=list(csv.DictReader(f))
    ntc={r['guide'] for r in guides if r['source_control_type']=='non_targeting_control'}
    if len(guides)!=1196 or len(ntc)!=106:raise ValueError('Wrong guide/control inventory')
    identity=frozen_identity(task,candidates);identity_sha=canonical_hash(identity)
    output.mkdir(parents=True,exist_ok=True);cache=output/'cache';cache.mkdir(exist_ok=True)
    checkpoint=output/'cache_checkpoint.json'
    saved=json.loads(checkpoint.read_text()) if checkpoint.exists() else dict(identity=identity,regions=[])
    if saved['identity']!=identity:raise ValueError('Checkpoint does not match exact input contracts')
    regions=candidates['regions']
    expected_regions={r['region_id']:r for r in regions}
    receipts={}
    for item in saved['regions']:
        rid=item['region']['region_id']
        if rid not in expected_regions or rid in receipts:raise ValueError('Unexpected/duplicate cached region')
        path=cache/(rid+'.jsonl.gz')
        validate_region_receipt(item,identity,expected_regions[rid],path);receipts[rid]=item
    total_bytes=sum(x['bytes'] for x in receipts.values())
    with pysam.AlignmentFile(task['source_BAM'],'rb',index_filename=task['source_index'],threads=1) as bam:
        ordered=[dict(name=n,length=int(l)) for n,l in zip(bam.references,bam.lengths)]
        if ordered!=task['reference_contigs_ordered']:raise ValueError('Full ordered reference inventory differs')
        if hashlib.sha256(bam.text.encode()).hexdigest()!=task['source_BAM_raw_header_sha256']:
            raise ValueError('Raw BAM header identity differs')
        with pysam.FastaFile(task['reference_fasta']) as ref:
            if [dict(name=n,length=int(l)) for n,l in zip(ref.references,ref.lengths)]!=ordered:
                raise ValueError('Full reference FASTA inventory differs')
        header=output/'source_header.sam'
        if header.exists() and header.read_text()!=bam.text:raise ValueError('Saved source header differs')
        if not header.exists():header.write_text(bam.text)
        for number,region in enumerate(regions):
            rid=region['region_id']
            if rid in receipts:continue
            final=cache/(rid+'.jsonl.gz');partial=cache/(rid+'.partial.gz')
            if final.exists() or partial.exists():
                raise ValueError('Unaccepted region artifact exists; preserve and diagnose before recovery')
            start=time.monotonic();fetched=selected=missing_umi=0;capped=False
            progress('extract_region',region=rid,region_index=number,total_regions=len(regions),
                     completed_regions=len(receipts),completed_guides=sum(len(x['region']['guides']) for x in receipts.values()))
            with partial.open('wb') as raw:
                with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as gz:
                    def emit(value):gz.write((json.dumps(value,separators=(',',':'))+'\n').encode())
                    emit(dict(kind='RNA_SAM_cache_header',identity=identity,region=region))
                    for read in bam.fetch(region['chrom'],region['start0'],region['end0']):
                        fetched+=1
                        if fetched>candidates['policy']['max_fetched_records_per_region']:
                            capped=True;break
                        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.mapping_quality<30:
                            continue
                        sample=samples_for_read(read,region,assigned,ntc)
                        if sample is None:continue
                        selected+=1;missing_umi+=int(not read.has_tag('UB'))
                        emit(dict(kind='original_SAM_record',sample=sample,sam=read.to_string()))
                raw.flush();os.fsync(raw.fileno())
            nbytes=partial.stat().st_size
            if total_bytes+nbytes>candidates['policy']['max_cache_bytes']:
                raise ValueError('Task cache safety bound reached; preserve completed regions and partial')
            partial.replace(final)
            item=dict(identity_sha256=identity_sha,region=region,
                      file='cache/'+final.name,bytes=nbytes,sha256=sha_file(final),
                      fetched_records=fetched,selected_records=selected,selected_without_UB=missing_umi,
                      scan_capped=capped,seconds=time.monotonic()-start)
            receipts[rid]=item;total_bytes+=nbytes
            atomic(checkpoint,dict(identity=identity,regions=[receipts[r['region_id']] for r in regions if r['region_id'] in receipts]))
            progress('region_cached',region=rid,completed_regions=len(receipts),total_regions=len(regions),
                     completed_guides=sum(len(x['region']['guides']) for x in receipts.values()),
                     total_guides=777,cache_bytes=total_bytes,selected_records=selected,scan_capped=capped)
    check_sources(task)
    assert len(receipts)==len(regions)
    summary=dict(task_id=task['task_id'],producer_commit=task['code_version'],identity_sha256=identity_sha,
                 cached_regions=len(receipts),cached_guides=777,cache_bytes=total_bytes,
                 selected_records=sum(x['selected_records'] for x in receipts.values()),
                 capped_regions=sum(x['scan_capped'] for x in receipts.values()),
                 stage='cache_complete_pending_allele_analysis',variants_analyzed=False,
                 DNA_or_functional_KO_efficiency_assessed=False)
    atomic(output/'cache_summary.json',summary)
    progress('cache_complete_pending_allele_analysis',cached_regions=len(receipts),
             cached_guides=777,cache_bytes=total_bytes,capped_regions=summary['capped_regions'])
    return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    config=json.loads(args.config.read_text());task=json.loads(Path(config['task_manifest']).read_text())
    if os.getuid()!=task['expected_owner_uid'] or socket.gethostname()!=task['verified_actual_hostname']:
        raise ValueError('Execution owner or host differs before any output write')
    if task['canonical_output_root']!='/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261030_all_target_editing_r1':
        raise ValueError('Output namespace differs before any write')
    candidates_path=Path(config['candidate_contract'])
    if sha_file(candidates_path)!=task['candidate_contract_sha256']:raise ValueError('Candidate manifest differs')
    candidates=json.loads(candidates_path.read_text())
    output=Path(task['canonical_output_root']);output.mkdir(parents=True,exist_ok=True)
    lock=output/'active_cache.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    s=Path('/proc/self/stat').read_text();startticks=int(s[s.rfind(')')+2:].split()[19])
    os.write(fd,json.dumps(dict(task_id=task['task_id'],pid=os.getpid(),uid=os.getuid(),
                               startticks=startticks,code_commit=task['code_version'],argv=sys.argv)).encode());os.close(fd)
    def progress(stage,**kw):
        value=dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),task_id=task['task_id'],**kw)
        atomic(output/'cache_status.json',value);print(json.dumps(value),flush=True)
    try:
        run(task,candidates,config['assignments'],config['guides'],progress)
    except BaseException as exc:
        progress('cache_failed_pending_recovery',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':main()
