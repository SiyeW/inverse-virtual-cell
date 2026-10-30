"""Extract a count-eligible PC checkpoint without reading expression matrix X."""
import argparse,hashlib,json,os,socket
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from rna_state_overview import select_records

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def prepare(cfg):
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    dest=Path(cfg['output']);assert not dest.exists()
    before={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in cfg['input_sha256']}
    for p,want in cfg['input_sha256'].items():
        if sha(p)!=want:raise ValueError('Accepted checkpoint or metadata differs: '+p)
    cells=pd.read_csv(cfg['cells'],dtype={'cell_barcode':str,'guide_id':str})
    assert len(cells)==cfg['expected_source_records']
    with h5py.File(cfg['checkpoint'],'r') as f:
        obs=read_elem(f['obs']);var=read_elem(f['var'])
        features=pd.read_csv(cfg['features'],index_col=0)
        assert list(features.index.astype(str))==list(var.index.astype(str))
        selected,rows=select_records(cells,obs.index.astype(str).tolist(),10)
        aligned=obs.loc[selected.cell_barcode]
        for col in ['guide_id','role']:
            assert aligned[col].astype(str).tolist()==selected[col].astype(str).tolist()
        assert np.array_equal(aligned.total_counts.to_numpy(),selected.total_counts.to_numpy())
        pc=f['obsm/X_pca'][()]
        assert pc.shape==(cfg['expected_checkpoint_records'],40) and np.isfinite(pc).all()
        arrays=dict(PCA=pc[rows].astype(np.float32),cell_barcode=selected.cell_barcode.to_numpy(dtype=str),guide=selected.guide_id.to_numpy(dtype=str),
                    role=selected.role.to_numpy(dtype=str),total_RNA=selected.total_counts.to_numpy(dtype=np.int64),RNA_ids=var.index.to_numpy(dtype=str))
    assert len(selected)==cfg['expected_selected_records']
    assert before=={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in before}
    temporary=dest.with_suffix('.receiving')
    with temporary.open('wb') as f:np.savez_compressed(f,**arrays)
    temporary.replace(dest)
    proof=dict(input_scope=cfg['input_scope'],output=str(dest),sha256=sha(dest),bytes=dest.stat().st_size,source_identity_before_after=before,
        records=len(selected),guides=selected.guide_id.nunique(),NTC_records=int(selected.role.eq('non_targeting_control').sum()),
        NTC_guides=selected.loc[selected.role.eq('non_targeting_control'),'guide_id'].nunique(),
        RNA_inventory_sha256=hashlib.sha256(json.dumps(var.index.astype(str).tolist()).encode()).hexdigest(),
        normalized_X_payload_read=False,original_counts_or_BAM_read=False,source_modified=False)
    Path(cfg['proof']).write_text(json.dumps(proof,indent=2))
    print(json.dumps(proof))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();prepare(json.loads(a.config.read_text()))
