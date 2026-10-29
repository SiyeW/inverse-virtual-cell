"""Accept a complete, source-verified cross-host RNA cache; never rename it BAM."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,socket,zipfile
from cache_guided_rna_alignments import atomic,canonical_hash,sha_file,frozen_identity


def validate_manifest_inventory(manifest,regions,identity):
    receipts=manifest['regions']
    expected={'source_header.sam'}|{'cache/'+rid+'.jsonl.gz' for rid in regions}
    if set(manifest['files'])!=expected:raise ValueError('Imported file roles differ')
    if len(receipts)!=len(regions) or {r['region']['region_id'] for r in receipts}!=set(regions):raise ValueError('Incomplete region inventory')
    for item in receipts:
        rid=item['region']['region_id']
        if item['region']!=regions[rid] or item['identity_sha256']!=canonical_hash(identity):raise ValueError('Region source reference inventory differs')
        if item['file']!='cache/'+rid+'.jsonl.gz':raise ValueError('Region file role differs')
        proof=manifest['files'][item['file']]
        if any(item[k]!=proof[k] for k in ['bytes','sha256']):raise ValueError('Receipt content inventory differs')
    return receipts


def main():
    p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True);p.add_argument('--candidates',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--archive',type=Path,required=True);a=p.parse_args()
    task=json.loads(a.task.read_text());candidate=json.loads(a.candidates.read_text());manifest=json.loads(a.manifest.read_text())
    if os.getuid()!=task['expected_owner_uid'] or socket.gethostname()!=task['verified_actual_hostname']:raise ValueError('Import execution identity differs')
    if task.get('dataset_id')!='chip01' or task.get('source_input_mode')!='verified_imported_RNA_cache_from_Windows_H':raise ValueError('Wrong imported artifact role')
    root=Path(task['canonical_output_root'])
    if str(root)!='/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261030_all_target_editing_r1/chip01':raise ValueError('Wrong source-specific canonical namespace')
    if sha_file(a.task)!=manifest['source_task_manifest_sha256'] or sha_file(a.candidates)!=task['candidate_contract_sha256']:raise ValueError('Source task or candidate identity differs')
    identity=frozen_identity(task,candidate)
    if manifest['identity']!=identity or not manifest['source_final_identity_verified']:raise ValueError('Incomplete or different source inventory')
    if manifest['source_host']!=task['source_host'] or manifest['source_owner_sid']!=task['source_owner_sid']:raise ValueError('Source host/owner differs')
    if sha_file(a.archive)!=manifest['archive_sha256']:raise ValueError('Imported archive differs')
    regions={r['region_id']:r for r in candidate['regions']};receipts=validate_manifest_inventory(manifest,regions,identity)
    root.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.archive) as z:
        if set(z.namelist())!=set(manifest['files']):raise ValueError('Archive member inventory differs')
        for name,item in manifest['files'].items():
            parts=Path(name).parts
            if Path(name).is_absolute() or '..' in parts:raise ValueError('Unsafe member path')
            b=z.read(name)
            if len(b)!=item['bytes'] or hashlib.sha256(b).hexdigest()!=item['sha256']:raise ValueError('Member content differs')
            target=root/name
            if target.exists():
                if sha_file(target)!=item['sha256']:raise ValueError('Different existing artifact remains unadopted')
            else:
                target.parent.mkdir(parents=True,exist_ok=True);temp=target.with_name(target.name+'.import_receiving');temp.write_bytes(b);temp.replace(target)
    if sha_file(root/'source_header.sam')!=task['source_BAM_raw_header_sha256']:raise ValueError('Raw source header differs')
    for item in receipts:
        rid=item['region']['region_id'];path=root/item['file']
        if item['identity_sha256']!=canonical_hash(identity) or item['region']!=regions[rid]:raise ValueError('Region source reference inventory differs')
        if path.stat().st_size!=item['bytes'] or sha_file(path)!=item['sha256']:raise ValueError('Region transport differs')
        with gzip.open(path,'rt') as f:header=json.loads(next(f))
        if header!=dict(kind='RNA_SAM_cache_header',identity=identity,region=regions[rid]):raise ValueError('Region biological identity is not equivalent')
    atomic(root/'cache_checkpoint.json',dict(identity=identity,regions=receipts))
    atomic(root/'import_receipt.json',dict(task_id=task['task_id'],dataset_id='chip01',
        import_manifest_sha256=sha_file(a.manifest),source_task_manifest_sha256=sha_file(a.task),
        identity_sha256=canonical_hash(identity),source_final_identity_verified=True,
        reference_contigs_ordered=task['reference_contigs_ordered'],source_gene_reference=task['source_gene_reference'],
        source_owner_sid=task['source_owner_sid'],source_host=task['source_host'],
        archive_sha256=manifest['archive_sha256'],regions=len(receipts)))
    atomic(root/'cache_summary.json',dict(task_id=task['task_id'],dataset_id='chip01',cached_regions=len(receipts),
        cached_guides=sum(len(r['region']['guides']) for r in receipts),cache_bytes=sum(r['bytes'] for r in receipts),
        capped_regions=sum(r['scan_capped'] for r in receipts),stage='cache_complete_pending_allele_analysis'))
    (Path(task['native_directory'])/'cache.exit').write_text('0\n')
    print(json.dumps(dict(imported=True,dataset_id='chip01',regions=len(receipts),
                         derived_cache_not_raw_BAM=True,complete_input_inventory_verified=True)))


if __name__=='__main__':main()
