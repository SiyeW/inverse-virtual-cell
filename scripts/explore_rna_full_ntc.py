"""Descriptive RNA-state exploration using complete query and NTC populations.

Energy U is the usual unequal-sample unbiased estimator. Its conditional
random-label variance is evaluated from pair-distance graph moments, with
tests against exhaustive permutations. This is not a guide-selection test.
All distances are streamed; no large distance matrix is saved or allocated.
"""
from pathlib import Path
import argparse,json,os,socket,sys,time
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances,pairwise_distances_chunked
from guide_rna_energy import atomic,sha,digest

LD=np.longdouble


def energy_u(n,m,within_query,cross,within_reference):
    if min(n,m)<2:raise ValueError('Both populations need at least two observations')
    return float(2*LD(cross)/(n*m)-LD(within_query)/(n*(n-1))-LD(within_reference)/(m*(m-1)))


def permutation_variance(N,n,total,ordered_square_sum,row_sum_square_sum):
    """Exact variance over every allocation of n of N fixed observations."""
    m=N-n
    if min(n,m)<2:raise ValueError('Invalid group sizes')
    def fall(x,k):return LD(0) if x<k else np.prod(np.asarray([x-j for j in range(k)],dtype=LD))
    a=-LD(2)/(n*(n-1));b=-LD(2)/(m*(m-1));c=LD(2)/(n*m)
    e2=(fall(n,2)*a*a+fall(m,2)*b*b+2*LD(n)*m*c*c)/fall(N,2)
    shared=(fall(n,3)*a*a+fall(m,3)*b*b+2*fall(n,2)*m*a*c+2*fall(m,2)*n*b*c+(LD(n)*fall(m,2)+LD(m)*fall(n,2))*c*c)/fall(N,3)
    disjoint=(fall(n,4)*a*a+fall(m,4)*b*b+2*fall(n,2)*fall(m,2)*a*b+4*fall(n,3)*m*a*c+4*fall(m,3)*n*b*c+4*fall(n,2)*fall(m,2)*c*c)/fall(N,4)
    edges=LD(total)/2;squares=LD(ordered_square_sum)/2
    overlap=LD(row_sum_square_sum)/2-squares;apart=(edges*edges-squares)/2-overlap
    variance=e2*squares+2*shared*overlap+2*disjoint*apart
    if variance < -LD('1e-8'):raise ArithmeticError('Negative permutation variance beyond rounding')
    return float(max(variance,0))


def moments(D):
    rows=D.sum(axis=1,dtype=np.float64)
    return rows.sum(dtype=LD),np.sum(D*D,dtype=LD),np.sum(rows.astype(LD)**2,dtype=LD)


class FullControl:
    def __init__(self,pc,working_memory=32,progress=None):
        self.pc=np.asarray(pc,dtype=np.float64);self.N=len(pc);self.memory=working_memory;self.row_sums=np.empty(self.N);self.square_sum=LD(0);start=0
        for D in pairwise_distances_chunked(self.pc,working_memory=working_memory,n_jobs=1):
            k=len(D);D[np.arange(k),start+np.arange(k)]=0
            self.row_sums[start:start+k]=D.sum(axis=1,dtype=np.float64);self.square_sum+=np.sum(D*D,dtype=LD);start+=k
            if progress and (start==self.N or start%2000<k):progress('NTC_distance_sums',completed=start,total=self.N)
        self.total=self.row_sums.sum(dtype=LD);self.row_square_sum=np.sum(self.row_sums.astype(LD)**2,dtype=LD)

    def evaluate(self,query,own_indices=None):
        q=np.asarray(query,dtype=np.float64);n=len(q);D=pairwise_distances(q);np.fill_diagonal(D,0)
        sx=D.sum(dtype=LD)
        if own_indices is not None:
            own=np.asarray(own_indices,dtype=int)
            if len(own)!=n or not np.array_equal(q,self.pc[own]):raise ValueError('Whole own control guide must be excluded')
            cross=self.row_sums[own].sum(dtype=LD)-sx;sy=self.total-2*self.row_sums[own].sum(dtype=LD)+sx;m=self.N-n
            variance=permutation_variance(self.N,n,self.total,self.square_sum,self.row_square_sum)
        else:
            m=self.N;cross=LD(0);cross_sq=LD(0);columns=np.zeros(m);query_rows=np.empty(n);start=0
            for C in pairwise_distances_chunked(q,self.pc,working_memory=self.memory,n_jobs=1):
                k=len(C);query_rows[start:start+k]=C.sum(axis=1,dtype=np.float64);columns+=C.sum(axis=0,dtype=np.float64);cross_sq+=np.sum(C*C,dtype=LD);start+=k
            cross=query_rows.sum(dtype=LD);sy=self.total
            combined_rows=np.concatenate([self.row_sums+columns,query_rows+D.sum(axis=1)])
            variance=permutation_variance(n+m,n,self.total+2*cross+sx,self.square_sum+2*cross_sq+np.sum(D*D,dtype=LD),np.sum(combined_rows.astype(LD)**2,dtype=LD))
        effect=energy_u(n,m,sx,cross,sy)
        return dict(query_records=n,reference_records=m,energy_unbiased=effect,random_label_SD=float(np.sqrt(variance)),
            within_query_sum=float(sx),cross_sum=float(cross),within_reference_sum=float(sy))


def run(cfg):
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    op=Path(cfg['native']);root=Path(cfg['output']);root.mkdir(exist_ok=True)
    for p,h in cfg['inputs_SHA256'].items():assert sha(p)==h,p
    fd=os.open(op/'active.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);s=Path('/proc/self/stat').read_text();ticks=int(s[s.rfind(')')+2:].split()[19]);os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,commit=cfg['commit'])).encode());os.close(fd)
    def progress(stage,**kw):atomic(op/'status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**kw))
    try:
        with np.load(cfg['pc_input'],allow_pickle=False) as z:pc=z['PCA'];labels=z['guide'];roles=z['role'];bc=z['cell_barcode']
        assert pc.shape==(71664,40) and len(set(bc))==71664
        old=pd.read_csv(cfg['identity_table']).set_index('guide',verify_integrity=True);groups=sorted(np.unique(labels));assert len(groups)==785
        control_idx=np.flatnonzero(roles=='non_targeting_control');assert len(control_idx)==19342
        control=FullControl(pc[control_idx],cfg['working_memory_MiB'],progress)
        nc_map=np.full(len(pc),-1,dtype=int);nc_map[control_idx]=np.arange(len(control_idx));contract=digest(cfg);rows=[];dest=root/'groups';dest.mkdir(exist_ok=True)
        for i,g in enumerate(groups):
            idx=np.flatnonzero(labels==g);role=str(roles[idx[0]]);assert len(np.unique(roles[idx]))==1 and len(idx)==int(old.loc[g,'n_each']) and role==old.loc[g,'role']
            path=dest/f'{i:04d}.json';members=digest(bc[idx].tolist())
            if path.exists():
                row=json.loads(path.read_text());assert row['contract_sha256']==contract and row['membership_SHA256']==members and row['guide']==g
            else:
                result=control.evaluate(pc[idx],nc_map[idx] if role=='non_targeting_control' else None)
                row=dict(guide=str(g),role=role,**result,contract_sha256=contract,membership_SHA256=members);atomic(path,row)
            rows.append(row)
            if (i+1)%10==0 or i+1==len(groups):progress('all_reference_effects',completed=i+1,total=len(groups))
        frame=pd.DataFrame(rows);frame.to_csv(root/'RNA_state_all_records.csv',index=False)
        grid=[]
        for n in range(10,int(frame.query_records.max())+1):
            v=permutation_variance(control.N,n,control.total,control.square_sum,control.row_square_sum);grid.append(dict(records=n,random_NTC_SD=float(np.sqrt(v))))
        pd.DataFrame(grid).to_csv(root/'NTC_random_label_variation.csv',index=False)
        summaries={}
        for kind,mask in [('targeting',frame.role.isin(['targeting','positive_control'])),('NTC',frame.role.eq('non_targeting_control')),('intergenic',frame.role.eq('intergenic_control'))]:
            f=frame.loc[mask];summaries[kind]=dict(groups=len(f),records=int(f.query_records.sum()),effect_quantiles=f.energy_unbiased.quantile([.1,.25,.5,.75,.9]).tolist(),SD_quantiles=f.random_label_SD.quantile([.1,.5,.9]).tolist())
        summary=dict(commit=cfg['commit'],records=len(pc),groups=len(frame),NTC_records=control.N,query_cells_downsampled=False,reference_NTC_downsampled=False,
            whole_own_NTC_guide_excluded=True,group_weight_in_distribution='one vote per observed guide, not one biological replicate',
            expression_space='unchanged pre-response normalized/log RNA, HVG2000 and PCA40',inference='descriptive RNA-state exploration, not a selected-guide list or KO efficiency',
            random_label_SD_is_not_biological_CI=True,comparisons=summaries,input_source_scope='company_merged_records_diagnostic_not_capture_resolved',contract_sha256=contract)
        atomic(root/'summary.json',summary);atomic(root/'manifest.json',{n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in ['RNA_state_all_records.csv','NTC_random_label_variation.csv','summary.json']});progress('complete_pending_plot_and_interpretation',completed=len(frame),total=len(frame))
    except BaseException as e:progress('failed_pending_recovery',error=repr(e));raise
    finally:(op/'active.lock').unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()))
