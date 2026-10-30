"""Prepare official SCEPTRE inputs from a verified invertible RNA checkpoint."""
import argparse,gzip,hashlib,json,os,socket,time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse,io
from threadpoolctl import threadpool_limits
from rna_state_overview import select_records
from recover_frozen_rna_counts import recover_counts

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def prepare(cfg):
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
    root=Path(cfg['output']);root.mkdir(exist_ok=True)
    proof_path=root/'input_proof.json';assert not proof_path.exists()
    for p,want in cfg['input_sha256'].items():
        if sha(p)!=want:raise ValueError('Accepted source differs: '+p)
    a=ad.read_h5ad(cfg['checkpoint'])
    source_cells=pd.read_csv(cfg['cells'],dtype={'cell_barcode':str,'guide_id':str})
    cells,rows=select_records(source_cells,a.obs_names.astype(str).tolist(),10)
    assert len(cells)==71664 and cells.guide_id.nunique()==785
    for col in ['guide_id','role']:
        assert a.obs.iloc[rows][col].astype(str).tolist()==cells[col].astype(str).tolist()
    qc=pd.read_csv(cfg['QC_table'],dtype={'cell_barcode':str,'guide_id':str}).set_index('cell_barcode',verify_integrity=True).reindex(cells.cell_barcode)
    assert qc.guide_id.astype(str).tolist()==cells.guide_id.astype(str).tolist()
    assert np.array_equal(qc.total_counts.to_numpy(),cells.total_counts.to_numpy())
    if not np.isfinite(qc.n_genes_by_counts).all() or np.any(qc.n_genes_by_counts<=0):raise ValueError('Invalid full-RNA detection counts')
    genes=a.var.copy();assert genes.index.is_unique and len(genes)==2000
    x=a.X[rows].tocsr()
    counts,recovery=recover_counts(x,cells.total_counts.to_numpy())
    assert counts.dtype==np.int32 and counts.shape==(71664,2000)
    catalogue=pd.read_csv(cfg['catalogue'],keep_default_na=False)
    needed=['guide','target_gene_name','target_gene_id','source_control_type']
    if not set(needed)<=set(catalogue.columns):raise ValueError('Missing design identity columns')
    catalogue=catalogue[needed].copy();assert len(catalogue)==1196 and catalogue.guide.is_unique
    guide_index={g:i for i,g in enumerate(catalogue.guide)}
    assigned=qc.guide_umis.to_numpy()
    if not np.isfinite(assigned).all() or np.any(assigned<1) or np.any(assigned!=np.rint(assigned)):raise ValueError('Invalid frozen assigned-guide UMI counts')
    guide_matrix=sparse.coo_matrix((assigned.astype('int32'),([guide_index[g] for g in cells.guide_id],np.arange(len(cells)))),shape=(len(catalogue),len(cells))).tocsr()
    meta=cells[['cell_barcode','guide_id','target','role']].copy()
    meta['total_RNA_full']=cells.total_counts.to_numpy();meta['detected_RNA_full']=qc.n_genes_by_counts.to_numpy()
    meta['assigned_guide_UMIs']=assigned
    meta.to_csv(root/'cells.tsv',sep='\t',index=False)
    feature=genes.reset_index();feature=feature.rename(columns={feature.columns[0]:'response_id'})
    if 'gene_name' not in feature:feature['gene_name']=feature.response_id
    feature.to_csv(root/'features.tsv',sep='\t',index=False)
    catalogue.to_csv(root/'guide_design.tsv',sep='\t',index=False)
    for name,matrix in [('response_counts.mtx.gz',counts.T),('assigned_guide_counts.mtx.gz',guide_matrix)]:
        temp=root/(name+'.receiving')
        with threadpool_limits(limits=2),gzip.open(temp,'wb',compresslevel=3) as f:io.mmwrite(f,matrix,field='integer',symmetry='general')
        temp.replace(root/name)
    names=['cells.tsv','features.tsv','guide_design.tsv','response_counts.mtx.gz','assigned_guide_counts.mtx.gz']
    proof=dict(created_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),input_scope='company_merged_records_diagnostic_not_capture_resolved',records=71664,response_RNAs=2000,
        observed_guides=785,design_guides=1196,NTC_records=int(meta.role.eq('non_targeting_control').sum()),
        original_counts_or_BAM_read=False,recovery=recovery,guide_matrix='frozen assigned-guide UMI counts only, not a new caller or complete raw capture matrix',
        expression_covariates='authoritative full-RNA totals and detected-RNA counts; not recomputed from the HVG subset',
        no_target_RNA_decrease_filter=True,files={n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in names})
    proof_path.write_text(json.dumps(proof,ensure_ascii=False,indent=2));print(json.dumps(proof))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);args=p.parse_args();prepare(json.loads(args.config.read_text()))
