"""Unsupervised overview of frozen RNA PCs before response-based selection.

Uses every QC/single-guide record meeting the proposed guide-count rule.
Guide identities never enter UMAP fitting. No response or KO filter is run.
"""
from pathlib import Path
import argparse,hashlib,json,os,socket,time,sys
import numpy as np
import pandas as pd


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def atomic(path,value):
    p=Path(path);tmp=p.with_suffix(p.suffix+'.receiving');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2));tmp.replace(p)


def select_records(cells,pc_barcodes,minimum=10):
    if not cells.cell_barcode.is_unique or len(set(pc_barcodes))!=len(pc_barcodes):raise ValueError('Duplicate cell identity')
    counts=cells.groupby('guide_id').size()
    selected=cells[cells.guide_id.map(counts).ge(minimum)].copy()
    position={str(b):i for i,b in enumerate(pc_barcodes)}
    missing=[b for b in selected.cell_barcode if b not in position]
    if missing:raise ValueError('Eligible records missing from frozen PCA')
    rows=np.array([position[b] for b in selected.cell_barcode],dtype=np.int64)
    return selected,rows


def calculate(cfg):
    import h5py
    from anndata.io import read_elem
    from importlib.metadata import version
    task=json.loads(Path(cfg['task_manifest']).read_text())
    assert os.getuid()==task['expected_owner_uid'] and socket.gethostname()==task['host']
    root=Path(task['canonical_output_root']);root.mkdir(parents=True,exist_ok=True)
    if root.exists() and (root/'summary.json').exists():raise ValueError('Completed output is immutable')
    lock=root/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    stat=Path('/proc/self/stat').read_text();startticks=int(stat[stat.rfind(')')+2:].split()[19])
    os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=startticks,argv=sys.argv,task_id=task['task_id'],commit=cfg['commit'])).encode());os.close(fd)
    def status(stage,**kw):atomic(root/'status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kw))
    def identity(p):
        s=Path(p).stat();return dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns)
    def verify():
        for p,expected in cfg['source_files'].items():
            if identity(p)!=dict(bytes=expected['bytes'],mtime_ns=expected['mtime_ns']):raise ValueError('RNA checkpoint/source metadata identity differs')
            if expected.get('sha256') and sha(p)!=expected['sha256']:raise ValueError('Source metadata content differs')
        if sha(__file__)!=cfg['script_sha256']:raise ValueError('Fixed execution source differs')
    try:
        verify();status('read_existing_pre_response_PCA_only')
        cache=root/'overview_input.npz';receipt_path=root/'overview_input_receipt.json'
        if cache.exists() or receipt_path.exists():
            if not (cache.exists() and receipt_path.exists()):raise ValueError('Incomplete input checkpoint; preserve for recovery')
            receipt=json.loads(receipt_path.read_text())
            if receipt['sha256']!=sha(cache) or receipt['config_sha256']!=cfg['config_sha256']:raise ValueError('Input checkpoint contract differs')
            with np.load(cache,allow_pickle=False) as z:arrays={k:z[k] for k in z.files}
        else:
            cells=pd.read_csv(task['input_cells'],dtype={'cell_barcode':str,'guide_id':str})
            if len(cells)!=72107:raise ValueError('Unexpected source cohort')
            with h5py.File(task['input_checkpoint'],'r') as f:
                obs=read_elem(f['obs']);var=read_elem(f['var'])
                selected,rows=select_records(cells,obs.index.astype(str).tolist(),10)
                aligned=obs.loc[selected.cell_barcode]
                for col in ['guide_id','role','target']:
                    if not np.array_equal(aligned[col].astype(str).to_numpy(),selected[col].astype(str).to_numpy()):raise ValueError('Barcode/guide/role alignment differs')
                if not np.array_equal(aligned.total_counts.to_numpy(),selected.total_counts.to_numpy()):raise ValueError('RNA totals differ')
                features=pd.read_csv(task['input_features'],index_col=0)
                if list(var.index.astype(str))!=list(features.index.astype(str)):raise ValueError('Ordered intended RNA inventory differs')
                pca=f['obsm/X_pca'][()]
                if pca.shape!=(72104,40) or not np.isfinite(pca).all():raise ValueError('Unexpected frozen PCA payload')
                pc=pca[rows].astype(np.float32);ratios=f['uns/pca/variance_ratio'][()].astype(float)
            if len(selected)!=71664:raise ValueError('Count-only eligibility changed')
            arrays=dict(PCA=pc,cell_barcode=selected.cell_barcode.to_numpy(dtype=str),guide=selected.guide_id.to_numpy(dtype=str),role=selected.role.to_numpy(dtype=str),total_RNA=selected.total_counts.to_numpy(dtype=np.int64),PCA_variance_ratio=ratios,RNA_ids=var.index.to_numpy(dtype=str))
            with (root/'overview_input.npz.receiving').open('wb') as f:np.savez_compressed(f,**arrays)
            (root/'overview_input.npz.receiving').replace(cache)
            receipt=dict(task_id=task['task_id'],config_sha256=cfg['config_sha256'],sha256=sha(cache),records=len(selected),source_PC_records=72104,original_QC_records=72107,
                whole_checkpoint_SHA256_is_historical_only=True,source_PC_array_SHA256=hashlib.sha256(pca.tobytes()).hexdigest(),ordered_RNA_ids_SHA256=hashlib.sha256(json.dumps(var.index.astype(str).tolist()).encode()).hexdigest(),
                excluded_QC_records_below10=72107-len(selected),no_response_selection=True,checkpoint_X_payload_read=False,original_count_payload_read=False)
            atomic(receipt_path,receipt)
        status('fit_unsupervised_UMAP',records=len(arrays['cell_barcode']),PCs=40,neighbors=15,min_dist=.5)
        import umap
        model=umap.UMAP(n_neighbors=15,min_dist=.5,n_components=2,metric='euclidean',random_state=20261030,n_jobs=1)
        coordinates=model.fit_transform(arrays['PCA']).astype(np.float32)
        if coordinates.shape!=(71664,2) or not np.isfinite(coordinates).all():raise ValueError('Invalid UMAP result')
        target=root/'overview_coordinates.npz'
        with (root/'overview_coordinates.npz.receiving').open('wb') as f:np.savez_compressed(f,UMAP=coordinates,**arrays)
        (root/'overview_coordinates.npz.receiving').replace(target)
        verify()
        summary=dict(task_id=task['task_id'],commit=cfg['commit'],records=71664,guide_count=len(np.unique(arrays['guide'])),counts_by_role={str(k):int(v) for k,v in zip(*np.unique(arrays['role'],return_counts=True))},
            retained_by_QC_and_guide_count_only=True,response_labels_used_for_fit=False,target_RNA_filter=False,normalization='existing total RNA to10000, log1p, Seurat HVG2000, PCA40 before any perturbation-signature calculation',
            UMAP=dict(n_neighbors=15,min_dist=.5,metric='euclidean',seed=20261030,n_jobs=1),PCA_variance_ratio=arrays['PCA_variance_ratio'].tolist(),
            versions={n:version(n) for n in ['numpy','pandas','anndata','h5py','umap-learn','numba','scikit-learn']},input_receipt=receipt,
            original_count_payload_read=False,checkpoint_X_payload_read=False,supervised_LDA_used=False,no_effect_threshold_selected=True)
        atomic(root/'summary.json',summary)
        files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [target,root/'summary.json',receipt_path]};atomic(root/'result_manifest.json',files)
        status('overview_complete_pending_plot_interpretation',records=71664)
    except BaseException as e:status('failed_pending_recovery',error=repr(e));raise
    finally:lock.unlink()


def plot(input_npz,summary_json,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','Noto Sans CJK SC','DejaVu Sans'],'axes.unicode_minus':False,'font.size':12,'svg.fonttype':'none'})
    with np.load(input_npz,allow_pickle=False) as z:xy=z['UMAP'];roles=z['role'];totals=z['total_RNA']
    summary=json.loads(Path(summary_json).read_text());assert len(xy)==summary['records']==71664
    colors=np.full(len(roles),'#ca9365',dtype='U7');colors[roles=='non_targeting_control']='#447d9e';colors[roles=='intergenic_control']='#939ba1'
    order=np.random.default_rng(20261030).permutation(len(xy))
    fig,ax=plt.subplots(1,2,figsize=(11.8,6.0),sharex=True,sharey=True)
    for a in ax:a.set_aspect('equal',adjustable='box');a.set_xlabel('UMAP 1');a.set_ylabel('UMAP 2');a.spines[['top','right']].set_visible(False)
    ax[0].scatter(xy[order,0],xy[order,1],c=colors[order],s=2,alpha=.4,linewidths=0,rasterized=True)
    labels=[('靶向sgRNA：52,128','#ca9365'),('NTC：19,342','#447d9e'),('基因间区对照：194','#939ba1')]
    ax[0].legend(handles=[Line2D([],[],marker='o',ls='',color=c,markersize=6,label=l) for l,c in labels],loc='upper right',frameon=False,fontsize=10)
    ax[0].set_title('NTC与靶向组分布在哪里？',pad=14)
    hue=np.log10(totals.astype(float));low,high=float(hue.min()),float(hue.max())
    a=ax[1].scatter(xy[order,0],xy[order,1],c=hue[order],s=2,cmap='viridis',vmin=low,vmax=high,alpha=.7,linewidths=0,rasterized=True)
    bar=fig.colorbar(a,ax=ax[1],fraction=.04,pad=.025);ticks=[np.log10(v) for v in [1000,2000,5000,10000,20000,50000,100000] if low<=np.log10(v)<=high];bar.set_ticks(ticks,[f'{round(10**v):,}' for v in ticks]);bar.set_label('每条记录测到的RNA总量（UMI）',fontsize=10)
    ax[1].set_title('同一张图：测到的RNA总量',pad=14)
    fig.suptitle('公司合并数据：筛选RNA响应前的整体RNA状态',fontsize=16,y=.98)
    fig.subplots_adjust(left=.07,right=.96,bottom=.2,top=.86,wspace=.24)
    fig.text(.5,.075,'每点一条细胞条码记录；只做既有RNA质控、单sgRNA分配及每sgRNA≥10条记录的数量筛选，共71,664条。\n两图坐标完全相同；未按靶RNA下降或RNA响应挑选记录，投影未使用sgRNA标签。图上重叠不等于统计上没有差异。',ha='center',va='center',fontsize=10,linespacing=1.7)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(out/('整体RNA状态概览.'+ext),dpi=300,facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('calculate');a.add_argument('--config',type=Path,required=True)
    b=sub.add_parser('plot');b.add_argument('--coordinates',type=Path,required=True);b.add_argument('--summary',type=Path,required=True);b.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.action=='calculate':calculate(json.loads(a.config.read_text()))
    else:plot(a.coordinates,a.summary,a.output)
