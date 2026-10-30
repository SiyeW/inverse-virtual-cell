"""Full-record unbiased energy estimates, without new RNA-count preprocessing."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

def estimates(pc,guide,role):
    control=role=='non_targeting_control';ids=np.flatnonzero(control);nt=pc[ids]
    d=cdist(nt,nt);rows=d.sum(axis=1);total=float(rows.sum());all_rows=[]
    for g in sorted(np.unique(guide)):
        mask=guide==g;n=int(mask.sum());r=np.unique(role[mask]);assert len(r)==1 and n>=10
        if r[0]=='non_targeting_control':
            q=np.flatnonzero(guide[ids]==g);within=float(d[np.ix_(q,q)].sum());cross=float(rows[q].sum()-within)
            other=total-2*cross-within;m=len(nt)-n
        else:
            x=pc[mask];within=float(cdist(x,x).sum());cross=float(cdist(x,nt).sum());other=total;m=len(nt)
        u=2*cross/(n*m)-within/(n*(n-1))-other/(m*(m-1))
        all_rows.append(dict(guide=g,role=r[0],records=n,reference_records=m,energy_unbiased=u))
    return pd.DataFrame(all_rows)

def run(cfg):
    p=Path(cfg['pc_input']);assert hashlib.sha256(p.read_bytes()).hexdigest()==cfg['pc_sha256']
    with np.load(p,allow_pickle=False) as z:
        assert z['PCA'].shape==(12591,40) and int((z['role']=='non_targeting_control').sum())==3769
        result=estimates(z['PCA'],z['guide'],z['role'])
    old=pd.read_csv(cfg['NTC_table']);old=old.drop_duplicates('guide').set_index('guide')
    nt=result[result.role.eq('non_targeting_control')].set_index('guide')
    assert set(nt.index)==set(old.index) and len(nt)==97
    error=float(np.max(np.abs(nt.energy_unbiased-old.loc[nt.index,'energy_U'])))
    assert error<1e-10,error
    result.to_csv(cfg['output'],index=False)
    Path(cfg['proof']).write_text(json.dumps(dict(all_guides=len(result),NTC_guides=97,NTC_bridge_max_abs_error=error,
        pc_input_sha256=cfg['pc_sha256'],count_or_BAM_read=False,all_frozen_records=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()))
