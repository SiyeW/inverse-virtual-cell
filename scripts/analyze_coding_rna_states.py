"""RNA-state exploration conditional on accepted project-ORF RNA evidence.

The existing target-only PCA coordinates and all records are preserved. Energy
tests use the full frozen log-normalized RNA panel, not the displayed PC axes.
Only guides with both sequence evidence classes enter a conditional comparison.
Cross-guide Pearson profiles use ALL such guides and ALL pairings, without an
RNA-response admission threshold, random control subset, or selected examples.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

from check_site_rna_association import energy_from_indicators,within_guide_permutations,adjusted_bh
from label_coding_rna_records import MODES,membership,display_category


def atomic(path,value):
    p=Path(path);temp=p.with_suffix(p.suffix+".writing")
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding="utf-8");temp.replace(p)


def comparison_masks(rows,guides,mode):
    flags=np.array([membership(row,mode) for row in rows],dtype=bool).reshape(-1,2)
    positive,negative=flags[:,0],flags[:,1]
    paired=[g for g in np.unique(guides) if positive[guides==g].any() and negative[guides==g].any()]
    use=(positive|negative)&np.isin(guides,paired)
    return positive,negative,use,paired


def residualize_amount(x,total):
    amount=np.log1p(np.asarray(total,float));z=np.column_stack([np.ones(len(amount)),amount-amount.mean()])
    return x-z@np.linalg.lstsq(z,x,rcond=None)[0]


def profile_correlation(left,right,genes,target_left,target_right):
    keep=(genes!=target_left)&(genes!=target_right)&np.isfinite(left)&np.isfinite(right)
    a,b=np.asarray(left[keep],float),np.asarray(right[keep],float)
    if len(a)<2:return np.nan
    a=a-a.mean();b=b-b.mean();aa=np.linalg.norm(a);bb=np.linalg.norm(b)
    # Do not report correlations of rounding noise after an exact small-n fit.
    if aa<1e-12 or bb<1e-12:return np.nan
    return float(np.clip(a@b/(aa*bb),-1,1))


def one_comparison(x,genes,target,rows,guides,total,mode,dataset):
    positive,negative,use,paired=comparison_masks(rows,guides,mode)
    n_positive=int((positive&use).sum());n_reference=int((negative&use).sum())
    result=dict(dataset=dataset,target_gene_id=target,comparison=mode,
                RNA_records=len(rows),all_frameshift_evidence_records=int(positive.sum()),
                all_matched_reference_records=int(negative.sum()),
                positive_records=n_positive,reference_records=n_reference,comparison_records=int(use.sum()),
                comparison_guides=len(paired),guides_with_frameshift=int(len(np.unique(guides[positive]))),
                status="too_few_observed_records")
    profiles=[]
    if not paired:
        result["status"]="no_within_guide_comparison" if positive.sum()>=2 and negative.sum()>=2 else "too_few_observed_records"
        return result,profiles
    xx=x[use].astype(float);gg=guides[use];yy=positive[use].astype(int);tt=np.asarray(total)[use]
    residual=residualize_amount(xx,tt)
    if mode=="internal_frameshift":
        for guide in paired:
            p=(gg==guide)&(yy==1);n=(gg==guide)&(yy==0)
            profiles.append(dict(guide=str(guide),target_gene_id=target,n_positive=int(p.sum()),n_reference=int(n.sum()),
                                 delta=xx[p].mean(axis=0)-xx[n].mean(axis=0),
                                 delta_amount_adjusted=residual[p].mean(axis=0)-residual[n].mean(axis=0)))
    if n_positive<2 or n_reference<2:return result,profiles
    seed=int.from_bytes(hashlib.sha256((dataset+target+mode).encode()).digest()[:4],"little")
    permutations,mutable=within_guide_permutations(yy,gg,repeats=999,seed=seed)
    if not mutable:raise ValueError("Paired guides unexpectedly have no exchangeable labels")
    keep=genes!=target
    with threadpool_limits(limits=2):
        for name,matrix in [("raw",xx),("RNA_amount_adjusted",residual)]:
            distances=cdist(matrix[:,keep],matrix[:,keep])
            observed=float(energy_from_indicators(distances,yy)[0])
            # A bounded batch avoids holding a full 999 x large-record matrix.
            exceed=0
            for offset in range(0,len(permutations),64):
                values=energy_from_indicators(distances,permutations[offset:offset+64])
                exceed+=int(np.count_nonzero(values>=observed-1e-12))
            result["energy_"+name]=observed
            result["p_"+name]=(1+exceed)/(1+len(permutations))
    result.update(status="tested",log_RNA_difference=float(np.log1p(tt[yy==1]).mean()-np.log1p(tt[yy==0]).mean()))
    return result,profiles


def all_profile_pairs(profiles,genes,dataset):
    rows=[]
    for i,left in enumerate(profiles):
        for right in profiles[i+1:]:
            row=dict(dataset=dataset,guide_left=left["guide"],guide_right=right["guide"],
                     target_left=left["target_gene_id"],target_right=right["target_gene_id"],
                     same_target=left["target_gene_id"]==right["target_gene_id"],
                     min_evidence_class_records=min(left["n_positive"],left["n_reference"],right["n_positive"],right["n_reference"]))
            for name,field in [("raw","delta"),("RNA_amount_adjusted","delta_amount_adjusted")]:
                row["correlation_"+name]=profile_correlation(left[field],right[field],genes,left["target_gene_id"],right["target_gene_id"])
            rows.append(row)
    return rows


def run(cfg):
    import anndata as ad
    import socket
    assert os.getuid()==2601 and socket.gethostname()=="rhesusbase"
    if hasattr(os,"sched_setaffinity"):
        os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
    labels_path=Path(cfg['labels']);labels=json.loads(labels_path.read_text(encoding='utf-8'))
    if hashlib.sha256(labels_path.read_bytes()).hexdigest()!=cfg['labels_sha256']:raise ValueError('Labels differ from accepted packet')
    catalogue=pd.read_csv(cfg['catalogue'],keep_default_na=False).set_index('guide')
    targets=labels['target_ids'];summaries=[];all_tests=[];all_pairs=[];pca_receipts={}
    out=Path(cfg["output"]);out.mkdir(parents=True,exist_ok=False)
    owner=dict(pid=os.getpid(),uid=os.getuid(),start_ticks=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
               config_sha256=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest())
    lock=out/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,json.dumps(owner).encode());os.close(fd)
    atomic(out/'owner.json',owner)
    try:
        for dataset in cfg['datasets']:
            did=dataset['id'];checkpoint=Path(dataset['RNA_checkpoint']);before=checkpoint.stat()
            if {'bytes':before.st_size,'mtime_ns':before.st_mtime_ns}!=dataset['checkpoint_identity']:
                raise ValueError('Frozen RNA checkpoint changed before reading')
            a=ad.read_h5ad(checkpoint);after=checkpoint.stat()
            assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
            if a.shape[1]!=2000 or not a.obs_names.is_unique:raise ValueError('Frozen RNA axes changed')
            meta=a.obs.copy();meta['target_gene_id']=meta.guide_id.astype(str).map(catalogue.target_gene_id)
            if meta.target_gene_id.isna().any():raise ValueError('Unmapped guide identity')
            genes=a.var_names.to_numpy(dtype=str)
            lookup={(r['guide'],r['cell_barcode']):r for r in labels[did]['records']}
            RNA_keys=set(zip(meta.guide_id.astype(str),meta.index.astype(str)))
            evidence_in_RNA=sum(k in RNA_keys for k in lookup)
            profiles=[];dataset_tests=[];target_records=0
            for i,target in enumerate(targets):
                ids=np.flatnonzero(meta.target_gene_id.eq(target).to_numpy())
                cells=meta.iloc[ids];guides=cells.guide_id.astype(str).to_numpy();barcodes=cells.index.to_numpy(dtype=str)
                totals=cells.total_counts.to_numpy(dtype=float);target_records+=len(ids)
                pca=Path(cfg['PCA_root'])/'RNA'/did/(target+'.npz')
                old=pca.stat();raw=pca.read_bytes();new=pca.stat()
                if (old.st_size,old.st_mtime_ns)!=(new.st_size,new.st_mtime_ns):raise ValueError('PCA cache changed during read')
                pca_receipts[did+'/'+target]=hashlib.sha256(raw).hexdigest()
                with np.load(io.BytesIO(raw),allow_pickle=False) as z:
                    if not (np.array_equal(z['guide'],guides) and np.array_equal(z['cell_barcode'],barcodes)
                            and np.array_equal(z['total_RNA'],totals)):
                        raise ValueError('Existing PCA and current RNA record identities differ')
                rows=[lookup.get((str(g),str(b)),{}) for g,b in zip(guides,barcodes)]
                categories=[display_category(r) for r in rows]
                count=lambda name:sum(c==name for c in categories)
                summary=dict(dataset=did,target_gene_id=target,target=str(catalogue.loc[catalogue.target_gene_id.eq(target),'target_gene_name'].iloc[0]),
                             RNA_records=len(rows),sgRNAs=len(set(guides)),
                             frameshift_records=count('Internal frameshift RNA'),
                             site_reference_records=count('Reference RNA at frameshift sites'),
                             other_variant_records=count('Other variant RNA'),
                             insufficient_site_records=count('Insufficient site information'),
                             frameshift_two_UMI_records=sum(r.get('frameshift_distinct_UMI_lower_bound',0)>=2 for r in rows))
                summary['frameshift_guides']=len(set(guides[np.array([c=='Internal frameshift RNA' for c in categories])]))
                if sum(summary[k] for k in ['frameshift_records','site_reference_records','other_variant_records','insufficient_site_records'])!=len(ids):
                    raise ValueError('Plot categories do not cover all records')
                x=a.X[ids].toarray()
                for mode in MODES:
                    result,guide_profiles=one_comparison(x,genes,target,rows,guides,totals,mode,did)
                    result['target']=summary['target'];dataset_tests.append(result)
                    if mode=='internal_frameshift':
                        summary['comparison_guides']=result['comparison_guides']
                        profiles.extend(guide_profiles)
                summaries.append(summary)
                atomic(out/'status.json',dict(stage='RNA_state_checks',dataset=did,completed=i+1,total=len(targets),target=target))
                print(json.dumps(dict(dataset=did,target=target,completed=i+1,total=len(targets))),flush=True)
            frame=pd.DataFrame(dataset_tests)
            for mode in MODES:
                for metric in ['raw','RNA_amount_adjusted']:
                    pcol='p_'+metric
                    if pcol not in frame:frame[pcol]=np.nan
                    valid=frame.comparison.eq(mode)&frame[pcol].notna()
                    frame.loc[valid,'q_'+metric]=adjusted_bh(frame.loc[valid,pcol].to_numpy())
            frame.to_csv(out/(did+'_association.csv'),index=False)
            all_tests.extend(frame.to_dict('records'))
            pairs=all_profile_pairs(profiles,genes,did);all_pairs.extend(pairs)
            pd.DataFrame(pairs,columns=['dataset','guide_left','guide_right','target_left','target_right','same_target',
                                       'min_evidence_class_records','correlation_raw','correlation_RNA_amount_adjusted']).to_csv(out/(did+'_guide_pairs.csv'),index=False)
            profile_meta=[{k:v for k,v in p.items() if k not in ('delta','delta_amount_adjusted')} for p in profiles]
            atomic(out/(did+'_profile_metadata.json'),profile_meta)
            if profiles:
                np.savez_compressed(out/(did+'_guide_profiles.npz'),gene=genes,
                                    delta=np.stack([p['delta'] for p in profiles]),
                                    delta_amount_adjusted=np.stack([p['delta_amount_adjusted'] for p in profiles]))
            atomic(out/(did+'_RNA_source_receipt.json'),dict(path=str(checkpoint),identity=dataset['checkpoint_identity'],
                   historical_accepted_sha256=dataset['historical_checkpoint_sha256'],
                   source_identity_unchanged_during_required_read=True,original_count_container_read=False,
                   RNA_features=2000,target_records=target_records,NTC_records=int(meta.role.eq('non_targeting_control').sum()),
                   source_sequence_evidence_records=len(lookup),sequence_evidence_records_in_frozen_RNA=evidence_in_RNA,
                   sequence_evidence_records_outside_frozen_RNA=len(lookup)-evidence_in_RNA))
            del a
        pd.DataFrame(summaries).to_csv(out/'case_summary.csv',index=False)
        atomic(out/'PCA_input_receipt.json',dict(source=cfg['PCA_root'],existing_coordinates_reused=True,sha256_during_required_read=pca_receipts))
        atomic(out/'method.json',dict(code_commit=cfg['commit'],labels_sha256=cfg['labels_sha256'],
              representation='Existing 2000 log-normalized RNAs; target RNA removed for energy tests and both target RNAs removed for each guide-profile pair.',
              comparison='Internal frameshift RNA versus reference RNA at the same guide candidate sites. Only guides with BOTH classes enter; all records remain in PCA plots.',
              null='999 evidence-label permutations within sgRNA; full observed sample sizes retained.',
              significance='BH separately within dataset, sensitivity and RNA-amount representation; exploratory associations, no KO selection.',
              sensitivities=list(MODES[1:]),RNA_amount_sensitivity='Regress RNA values on intercept and log1p(authoritative full-RNA total), without using evidence labels.',
              guide_profiles='Within each guide, all frameshift-evidence records mean minus all site-reference records mean, both classes >=1. All same/different target guide pairings, no significance gate or random subset.',
              pair_correlations='Descriptive Pearson r, not independent replicate tests; unequal record counts and company capture ambiguity remain.',
              PCA_recomputed=False,NTC='All prior NTC projections retained; NTC are not reference-allele cells in the evidence comparison.'))
        atomic(out/'status.json',dict(stage='statistics_complete_pending_figures',datasets=2,targets_per_dataset=104,
                                     records=sum(r['RNA_records'] for r in summaries),case_count=len(summaries)))
    except BaseException as error:
        atomic(out/'status.json',dict(stage='pending_recovery',error=repr(error)));raise
    finally:
        lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True)
    run(json.loads(parser.parse_args().config.read_text(encoding='utf-8')))
