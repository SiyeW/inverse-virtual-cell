"""Whole-cohort RNA allele observations from task-bound SAM caches.

bcftools 1.19 normalizes observed alleles and supplies an orthogonal candidate
call. Its diploid GT is never interpreted as a DNA genotype or KO efficiency.
Counts are corrected-CB/UB/GX families; conflicting PCR observations remain
ambiguous. Splice N operations never become deletion alleles.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time
import os
import socket
from functools import lru_cache

from cache_guided_rna_alignments import atomic,canonical_hash,sha_file,validate_region_receipt,frozen_identity,check_sources


@dataclass
class Read:
    sam: str
    sample: str

    def __post_init__(self):
        a=self.sam.split('\t')
        if len(a)<11:raise ValueError('Invalid SAM record')
        self.fields=a;self.flag=int(a[1]);self.chrom=a[2];self.start=int(a[3])-1
        self.seq=a[9].upper();self.qualities=[ord(q)-33 for q in a[10]] if a[10]!='*' else []
        self.tags={x.split(':',2)[0]:x.split(':',2)[2] for x in a[11:]}
        self.ops=[(op,int(n)) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',a[5])]
        if ''.join(f'{n}{op}' for op,n in self.ops)!=a[5]:raise ValueError('Unsupported CIGAR')
        if sum(n for op,n in self.ops if op in 'MIS=X')!=len(self.seq):raise ValueError('CIGAR/query length differs')

    @property
    def family(self):
        if not self.tags.get('CB') or not self.tags.get('UB'):return None
        return (self.sample,self.tags['CB'],self.tags['UB'],self.tags.get('GX',''))

    @property
    def signature(self):
        return (self.family,self.start,self.fields[5],self.seq,self.flag&16)

    @property
    def quality(self):return sum(self.qualities)

    def cohort_sam(self):
        a=self.fields[:11]+[x for x in self.fields[11:] if not x.startswith('RG:')]
        a.append('RG:Z:'+self.sample)
        return '\t'.join(a)+'\n'


def ref_piece(fetch,chrom,start,end):
    if start<0:return None
    value=fetch(chrom,start,end).upper()
    return value if len(value)==end-start and not (set(value)-set('ACGT')) else None


def raw_events(read,region,fetch,min_bq=20):
    if not read.qualities:return []
    cores=[(g['start0'],g['end0']) for g in region['guides']]
    def guides_for(start,end,insertion=False):
        return [g['guide'] for g in region['guides'] if
                (g['start0']<=start<=g['end0'] if insertion else start<g['end0'] and end>g['start0'])]
    events=[];r=read.start;q=0
    for op,n in read.ops:
        if op in 'M=X':
            intervals=sorted(set((max(r,a),min(r+n,b)) for a,b in cores if r<b and r+n>a))
            visited=set()
            for a,b in intervals:
                ref=fetch(read.chrom,a,b).upper()
                for i,base in enumerate(ref):
                    pos=a+i
                    if pos in visited:continue
                    visited.add(pos);qi=q+pos-r;alt=read.seq[qi]
                    if base in 'ACGT' and alt in 'ACGT' and alt!=base and read.qualities[qi]>=min_bq:
                        events.append(dict(chrom=read.chrom,pos0=pos,ref=base,alt=alt,guides=guides_for(pos,pos+1),kind='SNV'))
            r+=n;q+=n
        elif op=='I':
            ids=guides_for(r,r,True)
            flank=read.qualities[max(0,q-1):min(len(read.seq),q+n+1)]
            anchor=ref_piece(fetch,read.chrom,r-1,r)
            alt=read.seq[q:q+n]
            if ids and anchor and len(flank)==n+2 and min(flank)>=min_bq and not (set(alt)-set('ACGT')):
                events.append(dict(chrom=read.chrom,pos0=r-1,ref=anchor,alt=anchor+alt,guides=ids,kind='INS'))
            q+=n
        elif op=='D':
            ids=guides_for(r,r+n);flank=read.qualities[max(0,q-1):min(len(read.seq),q+1)]
            ref=ref_piece(fetch,read.chrom,r-1,r+n) if n<=1000 else None
            if ids and ref and len(flank)==2 and min(flank)>=min_bq:
                events.append(dict(chrom=read.chrom,pos0=r-1,ref=ref,alt=ref[0],guides=ids,kind='DEL'))
            r+=n
        elif op=='N':r+=n
        elif op=='S':q+=n
    return events


def event_key(e):return (e['chrom'],int(e['pos0']),e['ref'],e['alt'])


def covers_reference(read,start,end,min_bq=20):
    """A fully aligned high-quality reference window; N/D/I aren't WT bases."""
    r=read.start;q=0;seen={}
    for op,n in read.ops:
        if op in 'M=X':
            lo=max(r,start);hi=min(r+n,end)
            if lo<hi:
                for pos in range(lo,hi):
                    qi=q+pos-r
                    if not read.qualities or read.qualities[qi]<min_bq:return None
                    seen[pos]=read.seq[qi]
            r+=n;q+=n
        elif op in 'DN':
            if r<end and r+n>start:return None
            r+=n
        elif op=='I':
            if start<r<=end:return None
            q+=n
        elif op=='S':q+=n
    return ''.join(seen[p] for p in range(start,end)) if len(seen)==end-start else None


def contextual_flags(read,start,end):
    r=read.start;skips=False;soft=False
    for op,n in read.ops:
        if op=='N' and r<end and r+n>start:skips=True
        if op=='S' and start-5<=r<=end+5:soft=True
        if op in 'MDN=X':r+=n
    return skips,soft


def has_observed_base(read,start,end,min_bq=20):
    """Equivalent to the old per-position covers_reference predicate, faster."""
    if not read.qualities:return False
    r=read.start;blocked=set()
    for op,n in read.ops:
        if op=='I':blocked.add(r-1)
        if op in 'MDN=X':r+=n
    r=read.start;q=0
    for op,n in read.ops:
        if op in 'M=X':
            a=max(start,r);b=min(end,r+n)
            if any(p not in blocked and read.qualities[q+p-r]>=min_bq for p in range(a,b)):return True
            r+=n;q+=n
        elif op in 'DN':r+=n
        elif op in 'IS':q+=n
    return False


def normalize_events(events,reference,header,bcftools,directory):
    source=directory/'observed_events.vcf';normal=directory/'normalized_events.vcf'
    lookup={};lines=['##fileformat=VCFv4.2']
    for x in header:lines.append(f'##contig=<ID={x["name"]},length={x["length"]}>')
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    order={x['name']:i for i,x in enumerate(header)}
    keys=sorted(events,key=lambda k:(order[k[0]],k[1],k[2],k[3]))
    for i,k in enumerate(keys):
        name=f'E{i:07d}';lookup[name]=k
        lines.append(f'{k[0]}\t{k[1]+1}\t{name}\t{k[2]}\t{k[3]}\t.\t.\t.')
    source.write_text('\n'.join(lines)+'\n')
    command=[bcftools,'norm','-f',reference,'-c','e','-Ov','-o',str(normal),str(source)]
    p=subprocess.run(command,capture_output=True,text=True,timeout=60)
    (directory/'norm.stderr').write_text(p.stderr)
    if p.returncode:raise ValueError('bcftools reference normalization failed: '+p.stderr[-1000:])
    result={}
    for line in normal.read_text().splitlines():
        if line.startswith('#'):continue
        a=line.split('\t');raw=lookup[a[2]]
        result[raw]=(a[0],int(a[1])-1,a[3],a[4])
    if len(result)!=len(keys):raise ValueError('Normalization lost candidate identities')
    return result


def bcftools_candidates(reads,region,task,reference,bcftools,directory):
    # One actual original read per corrected UMI family. Conflicts are handled
    # separately by family evidence, not resolved by the candidate genotype.
    representatives={}
    for read in reads:
        f=read.family
        if f is not None and (f not in representatives or read.quality>representatives[f].quality):representatives[f]=read
    sam=directory/'candidate_input.sam'
    samples=sorted({r.sample for r in representatives.values()})
    with sam.open('w') as f:
        f.write('@HD\tVN:1.6\tSO:coordinate\n')
        for x in task['reference_contigs_ordered']:f.write(f'@SQ\tSN:{x["name"]}\tLN:{x["length"]}\n')
        for s in samples:f.write(f'@RG\tID:{s}\tSM:{s}\n')
        f.write('@CO\tDerived RNA analysis cohorts, not biological DNA samples; diploid GT must not be interpreted.\n')
        for read in sorted(representatives.values(),key=lambda r:(r.start,r.fields[0])):f.write(read.cohort_sam())
    if not representatives:return {}
    mp=directory/'pileup.vcf';called=directory/'caller_candidates.vcf'
    commands=[
        [bcftools,'mpileup','-f',reference,'-q','30','-Q','20','-D','-d','100000','-L','100000','-m','2','-p',
         '--ns','UNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY','-a','FORMAT/AD,FORMAT/DP','-Ov','-o',str(mp),str(sam)],
        [bcftools,'call','-m','-v','-Ov','-o',str(called),str(mp)],
    ]
    logs=[]
    for command in commands:
        p=subprocess.run(command,capture_output=True,text=True,timeout=180)
        logs.append(dict(command=command,code=p.returncode,stderr=p.stderr))
        if p.returncode:
            atomic(directory/'caller_log.json',logs);raise ValueError('bcftools candidate stage failed')
    atomic(directory/'caller_log.json',logs)
    normalized=directory/'caller_candidates.normalized.vcf'
    command=[bcftools,'norm','-a','-m','-any','-f',reference,'-c','e','-Ov','-o',str(normalized),str(called)]
    p=subprocess.run(command,capture_output=True,text=True,timeout=60)
    if p.returncode:raise ValueError('Caller candidate normalization failed: '+p.stderr[-1000:])
    observations={}
    for line in normalized.read_text().splitlines():
        if line.startswith('#'):continue
        a=line.split('\t')
        for alt in a[4].split(','):
            if alt in ['*','<*>','.']:continue
            observations[(a[0],int(a[1])-1,a[3],alt)]=dict(quality=None if a[5]=='.' else float(a[5]),filter=a[6])
    return observations


def analyze_region(cache,receipt,identity,region,task,reference,output,bcftools,analysis_contract_sha256):
    validate_region_receipt(receipt,identity,region,cache)
    output.mkdir(parents=True,exist_ok=True)
    if list(output.iterdir()):raise ValueError('Unaccepted partial analysis exists; archive and diagnose before recovery')
    all_reads={};raw_records=missing_umi=0
    with gzip.open(cache,'rt') as f:
        next(f)
        for line in f:
            item=json.loads(line);raw_records+=1
            if item['kind']!='original_SAM_record':raise ValueError('Unknown cache record schema')
            read=Read(item['sam'],item['sample'])
            if read.sample not in ['NTC']+[g['sample_id'] for g in region['guides']]:
                raise ValueError('Cached analysis-cohort identity differs')
            if read.family is None:missing_umi+=1;continue
            key=read.signature
            if key not in all_reads or read.quality>all_reads[key].quality:all_reads[key]=read
    reads=list(all_reads.values())
    import pysam
    events={};read_events=[]
    with pysam.FastaFile(reference) as fasta:
        fetch=lru_cache(maxsize=4096)(fasta.fetch)
        for read in reads:
            got=raw_events(read,region,fetch)
            keys=set()
            for e in got:
                key=event_key(e);keys.add(key)
                if key not in events:events[key]=dict(e,guides=set(e['guides']))
                else:events[key]['guides'].update(e['guides'])
            read_events.append(keys)
        normalized=normalize_events(events,reference,task['reference_contigs_ordered'],bcftools,output)
        norm_guides=defaultdict(set)
        for raw,key in normalized.items():norm_guides[key].update(events[raw]['guides'])
        by_family=defaultdict(list)
        for read,raw in zip(reads,read_events):by_family[read.family].append((read,{normalized[k] for k in raw}))
        variant_rows=[]
        for key,guide_ids in sorted(norm_guides.items()):
            chrom,pos,ref,alt=key
            sample_counts=defaultdict(lambda:dict(alt=0,ref=0,other=0,ambiguous=0,alt_barcodes=set(),ref_barcodes=set()))
            for family,observations in by_family.items():
                states=set()
                for read,event_set in observations:
                    if key in event_set:states.add('alt');continue
                    observed=covers_reference(read,pos,pos+len(ref))
                    if observed is not None:states.add('ref' if observed==ref else 'other')
                if not states:continue
                state=next(iter(states)) if len(states)==1 else 'ambiguous'
                c=sample_counts[family[0]];c[state]+=1
                if state in ['alt','ref']:c[state+'_barcodes'].add(family[1])
            for guide in sorted(guide_ids):
                sample=next(g['sample_id'] for g in region['guides'] if g['guide']==guide)
                own=sample_counts[sample];ntc=sample_counts['NTC']
                row=dict(guide=guide,chrom=chrom,pos0=pos,ref=ref,alt=alt,
                         kind='SNV' if len(ref)==len(alt)==1 else ('DEL' if len(ref)>len(alt) else 'INS'),
                         length_change=len(alt)-len(ref),scan_capped=receipt['scan_capped'],
                         function_consequence='not_assigned_correct_target_ORF_not_yet_verified',
                         not_DNA_genotype=True)
                for prefix,counts in [('guide',own),('NTC',ntc)]:
                    for k,v in counts.items():row[prefix+'_'+k]=len(v) if isinstance(v,set) else v
                variant_rows.append(row)
        # Whole core opportunities are descriptive RNA family observations.
        guide_rows=[]
        for guide in region['guides']:
            own_families=set();ntc_families=set();own_spliced=set();own_softclip=set()
            for family,observations in by_family.items():
                if family[0] not in [guide['sample_id'],'NTC']:continue
                present=False
                for read,_ in observations:
                    skipped,clipped=contextual_flags(read,guide['start0'],guide['end0'])
                    if family[0]==guide['sample_id']:
                        if skipped:own_spliced.add(family)
                        if clipped:own_softclip.add(family)
                    if has_observed_base(read,guide['start0'],guide['end0']):present=True
                    if present:break
                if present:(own_families if family[0]==guide['sample_id'] else ntc_families).add(family)
            owned=[r for r in variant_rows if r['guide']==guide['guide']]
            guide_rows.append(dict(guide=guide['guide'],region=region['region_id'],own_RNA_families_observed=len(own_families),
                                   NTC_RNA_families_observed=len(ntc_families),candidate_variants=len(owned),
                                   own_RNA_families_with_splice_skip=len(own_spliced),
                                   own_RNA_families_with_nearby_softclip=len(own_softclip),
                                   max_alt_families=max([r['guide_alt'] for r in owned],default=0),scan_capped=receipt['scan_capped']))
    called=bcftools_candidates(reads,region,task,reference,bcftools,output)
    for row in variant_rows:
        key=event_key(row);row['bcftools_candidate_seen']=key in called
        row['bcftools_candidate_quality']=called.get(key,{}).get('quality')
    result=dict(task_id=task['task_id'],region=region,source_cache_sha256=receipt['sha256'],
                frozen_source_identity_sha256=canonical_hash(identity),raw_cached_records=raw_records,
                analysis_contract_sha256=analysis_contract_sha256,
                missing_UMI_records=missing_umi,unique_read_observations=len(reads),RNA_families=len(by_family),
                variants=variant_rows,guides=guide_rows,bcftools_candidate_records=len(called),
                caller_GT_ignored=True,functional_KO_efficiency_assessed=False)
    atomic(output/'result.json',result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    cfg=json.loads(args.config.read_text())
    if sha_file(cfg['task_manifest'])!=cfg['producer_task_manifest_sha256']:
        raise ValueError('Producer task manifest bytes differ')
    if sha_file(__file__)!=cfg['analysis_script_sha256']:
        raise ValueError('Fixed analysis source bytes differ')
    task=json.loads(Path(cfg['task_manifest']).read_text())
    if os.getuid()!=task['expected_owner_uid'] or socket.gethostname()!=task['verified_actual_hostname']:
        raise ValueError('Wrong execution owner or host')
    allowed=['/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261030_all_target_editing_r1',
             '/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261030_all_target_editing_r1/chip01']
    if task['canonical_output_root'] not in allowed:
        raise ValueError('Wrong parent output namespace')
    def verify_inputs():
        if task.get('source_input_mode')=='verified_imported_RNA_cache_from_Windows_H':
            if task.get('dataset_id')!='chip01' or not cfg.get('import_receipt_sha256'):raise ValueError('Unapproved imported artifact role')
            receipt=Path(task['canonical_output_root'])/'import_receipt.json'
            if sha_file(receipt)!=cfg['import_receipt_sha256']:raise ValueError('Import receipt differs')
            imported=json.loads(receipt.read_text())
            if not imported['source_final_identity_verified'] or imported['source_task_manifest_sha256']!=cfg['producer_task_manifest_sha256']:
                raise ValueError('Imported cache lacks complete source verification')
            from cache_guided_rna_alignments import inventory
            if inventory(task['reference_fasta'])!=task['reference_identity'] or inventory(task['reference_index'])!=task['reference_index_identity']:
                raise ValueError('Analysis reference differs')
        else:check_sources(task)
    verify_inputs()
    root=Path(task['canonical_output_root']);out=root/'alleles';out.mkdir(exist_ok=True)
    # The consumer code is separately fixed; cache identity remains producer-bound.
    lock=out/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    stat=Path('/proc/self/stat').read_text();startticks=int(stat[stat.rfind(')')+2:].split()[19])
    os.write(fd,json.dumps(dict(task_id=task['task_id'],pid=os.getpid(),uid=os.getuid(),startticks=startticks,
                               analysis_commit=cfg['analysis_commit'],argv=list(os.sys.argv))).encode());os.close(fd)
    def progress(stage,**values):
        value=dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**values)
        atomic(out/'status.json',value);print(json.dumps(value),flush=True)
    try:
        version=subprocess.check_output([cfg['bcftools'],'--version'],text=True).splitlines()[0]
        if version!='bcftools 1.19':raise ValueError('Unexpected frozen caller version')
        if sha_file(cfg['bcftools'])!=cfg['bcftools_binary_sha256']:
            raise ValueError('Caller binary content differs')
        source=Path(cfg['candidate_contract'])
        if sha_file(source)!=task['candidate_contract_sha256']:raise ValueError('Candidate contract differs')
        candidates=json.loads(source.read_text());regions=candidates['regions'];done={}
        expected_identity=frozen_identity(task,candidates)
        analysis_contract=dict(task_id=task['task_id'],analysis_commit=cfg['analysis_commit'],
            producer_identity_sha256=canonical_hash(expected_identity),
            source_gene_reference=task['source_gene_reference'],bcftools_version='1.19',
            bcftools_binary_sha256=sha_file(cfg['bcftools']),
            method='RNA allele observations, corrected-family conflict exclusion, bcftools candidate-only GT ignored',
            quality=dict(min_MAPQ=30,min_BQ=20),norm_reference_sha256=task['reference_metadata_sha256'])
        contract_sha=canonical_hash(analysis_contract)
        old_contract=out/'analysis_contract.json'
        if old_contract.exists() and json.loads(old_contract.read_text())!=analysis_contract:
            raise ValueError('Previous analysis is not equivalent in source/gene/reference/code contracts')
        atomic(old_contract,analysis_contract)
        while True:
            checkpoint=root/'cache_checkpoint.json'
            if not checkpoint.exists():time.sleep(20);continue
            ledger=json.loads(checkpoint.read_text());identity=ledger['identity']
            if identity!=expected_identity:raise ValueError('Complete producer source/reference/assignment/code inventory differs')
            receipts={x['region']['region_id']:x for x in ledger['regions']}
            for region in regions:
                rid=region['region_id']
                if rid not in receipts or rid in done:continue
                result_path=out/rid/'result.json'
                if result_path.exists():
                    receipt_path=result_path.parent/'result_receipt.json'
                    if not receipt_path.exists():raise ValueError('Existing result lacks a completed artifact receipt')
                    result_receipt=json.loads(receipt_path.read_text())
                    if result_receipt['analysis_contract_sha256']!=contract_sha:
                        raise ValueError('Previous artifact contract differs')
                    for name,item in result_receipt['files'].items():
                        path=result_path.parent/name
                        if Path(name).name!=name or path.stat().st_size!=item['bytes'] or sha_file(path)!=item['sha256']:
                            raise ValueError('Saved analysis artifact bytes differ')
                    old=json.loads(result_path.read_text())
                    if old['source_cache_sha256']!=receipts[rid]['sha256'] or old['region']!=region:
                        raise ValueError('Existing allele artifact is not equivalent')
                    if old.get('analysis_commit')!=cfg['analysis_commit']:
                        raise ValueError('Existing result lacks exact analysis-code provenance')
                    if old.get('analysis_contract_sha256')!=contract_sha:
                        raise ValueError('Previous result lacks equivalent complete reference/gene/source contract')
                    result=old
                else:
                    progress('analyze_cached_region',region=rid,completed=len(done),total=len(regions))
                    result=analyze_region(root/receipts[rid]['file'],receipts[rid],identity,region,task,
                                          task['reference_fasta'],out/rid,cfg['bcftools'],contract_sha)
                    result['analysis_commit']=cfg['analysis_commit'];atomic(result_path,result)
                    files={p.name:dict(bytes=p.stat().st_size,sha256=sha_file(p)) for p in result_path.parent.iterdir() if p.is_file()}
                    atomic(result_path.parent/'result_receipt.json',dict(task_id=task['task_id'],analysis_contract_sha256=contract_sha,files=files))
                done[rid]=result
                progress('region_alleles_complete',region=rid,completed=len(done),total=len(regions))
            if len(done)==len(regions):break
            if (Path(cfg['producer_exit'])).exists():
                raise ValueError('Producer exited before complete region inventory; recovery required')
            time.sleep(20)
        guide_rows=[x for r in regions for x in done[r['region_id']]['guides']]
        variants=[x for r in regions for x in done[r['region_id']]['variants']]
        atomic(out/'all_guide_observations.json',guide_rows);atomic(out/'all_variant_observations.json',variants)
        verify_inputs()
        progress('allele_analysis_complete_pending_whole_cohort_interpretation',regions=len(done),guides=len(guide_rows),variants=len(variants))
    except BaseException as exc:progress('failed_pending_recovery',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':main()
