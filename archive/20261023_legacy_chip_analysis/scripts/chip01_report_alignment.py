"""Finite control-annotation and descriptive effect-size QC, not hit calling."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import sparse


TYPES = ['non_targeting_control', 'intergenic_control', 'positive_control', 'targeting']


def annotation_overlay(coverage, annotations):
    """Keep the frozen coarse labels while adding a fully covered source overlay."""
    a = pd.DataFrame(annotations)
    required = {'guide', 'current_target', 'current_type', 'comparator_type'}
    if not required.issubset(a) or a.guide.duplicated().any() or coverage.name.duplicated().any():
        raise ValueError('Unique complete catalogue annotations required')
    if set(a.guide.astype(str)) != set(coverage.name.astype(str)):
        raise ValueError('Catalogue overlay must cover the exact frozen guide universe')
    if not set(a.comparator_type).issubset(TYPES):
        raise ValueError('Unrecognized source control category')
    frozen = coverage.set_index('name').copy()
    a = a.set_index('guide').loc[frozen.index]
    for column, values in [('target_gene_name', a.current_target), ('guide_type', a.current_type)]:
        if column not in frozen or not np.array_equal(frozen[column].astype(str), values.astype(str)):
            raise ValueError('Overlay conflicts with frozen reference: '+column)
    expected_negative = a.comparator_type.isin(TYPES[:2]).to_numpy()
    if not np.array_equal(expected_negative, frozen.guide_type.eq('NTC').to_numpy()):
        raise ValueError('Coarse/source control identities conflict')
    frozen['source_control_type'] = a.comparator_type
    frozen['source_role_validated_in_this_assay'] = False
    return frozen.reset_index().rename(columns={'name':'guide'})


def reference_masks(obs, overlay, historical_core):
    labels = overlay.set_index('guide').source_control_type
    mapped = obs.guide_id.astype(str).map(labels)
    if mapped.isna().any():
        raise ValueError('Source cell has no unique catalogue annotation')
    historical_core = np.asarray(historical_core, dtype=bool)
    pure = mapped.eq('non_targeting_control').to_numpy()
    if historical_core.shape != pure.shape or np.any(historical_core & ~obs.is_control.to_numpy(dtype=bool)):
        raise ValueError('Historical reference/source identity differs')
    revised = historical_core & pure
    if revised.sum() < 2:
        raise ValueError('Insufficient pure NTC core reference')
    return mapped.to_numpy(dtype=str), pure, revised


def leave_guide_out_reference(normalized, guides, reference_mask, guide):
    """For a control guide, remove *all* its cells from its RNA reference."""
    reference = np.asarray(reference_mask, dtype=bool) & (np.asarray(guides,dtype=str) != guide)
    if reference.sum() < 2:
        raise ValueError('Insufficient nonoverlapping guide reference')
    mean = np.asarray(normalized[reference].mean(axis=0)).ravel()
    return mean, np.flatnonzero(reference)


def expression_effect(mean, reference, pseudocount=1.):
    mean, reference = np.asarray(mean,float), np.asarray(reference,float)
    if mean.shape != reference.shape or not np.isfinite(mean).all() or not np.isfinite(reference).all():
        raise ValueError('RNA mean/reference identity mismatch')
    if (mean<0).any() or (reference<0).any() or not np.isfinite(pseudocount) or pseudocount<=0:
        raise ValueError('Nonnegative RNA and positive declared pseudocount required')
    return np.log2((mean+pseudocount)/(reference+pseudocount))


def unique_direct_feature(target, selected_names, global_names):
    selected = np.flatnonzero(np.asarray(selected_names,dtype=str)==target)
    if len(selected)==1 and Counter(np.asarray(global_names,dtype=str))[target]==1:
        return int(selected[0])
    return None


def stable_rng(seed, key):
    token=int.from_bytes(sha256(str(key).encode()).digest()[:8],'big')
    return np.random.default_rng(np.random.SeedSequence([int(seed),token]))


def effect_profiles(normalized, obs, overlay, reference, feature_names, global_names,
                    min_cells=10, matched_cells=10, repeats=8, seed=0, pseudocount=1.):
    """All-cell and size-matched guide profiles; no p-values or cell filtering."""
    if min_cells<matched_cells or min_cells<1 or repeats<1:
        raise ValueError('Invalid fixed-size diagnostic budget')
    guides=obs.guide_id.astype(str).to_numpy()
    if normalized.shape!=(len(obs),len(feature_names)):
        raise ValueError('Selected RNA feature/row identity differs')
    means=[];reference_means=[];effects=[];matched_first=[];records=[];sampling=[];metric_rows=[]
    for row in overlay.sort_values('guide').itertuples():
        ix=np.flatnonzero(guides==row.guide)
        if len(ix)<min_cells:continue
        center,ref_ix=leave_guide_out_reference(normalized,guides,reference,row.guide)
        if np.intersect1d(ix,ref_ix).size:raise ValueError('Guide/reference overlap')
        mean=np.asarray(normalized[ix].mean(0)).ravel()
        profile=expression_effect(mean,center,pseudocount)
        direct=unique_direct_feature(str(row.target_gene_name),feature_names,global_names)
        if direct is not None:profile[direct]=np.nan
        records.append(dict(guide=row.guide,target=str(row.target_gene_name),cells=len(ix),
            source_control_type=row.source_control_type,reference_cells=len(ref_ix),
            direct_feature_index=direct,reference_overlap_cells=0))
        means.append(mean);reference_means.append(center);effects.append(profile)
        rng=stable_rng(seed,row.guide)
        for repeat in range(repeats):
            chosen=np.sort(rng.choice(ix,matched_cells,replace=False))
            sampled=np.asarray(normalized[chosen].mean(0)).ravel()
            effect=expression_effect(sampled,center,pseudocount)
            if direct is not None:effect[direct]=np.nan
            finite=effect[np.isfinite(effect)]
            metric_rows.append(dict(guide=row.guide,target=str(row.target_gene_name),
                source_control_type=row.source_control_type,repeat=repeat+1,cells=matched_cells,
                reference_cells=len(ref_ix),reference_overlap_cells=0,features=len(finite),
                rms_log2fc=float(np.sqrt(np.mean(finite**2))),median_log2fc=float(np.median(finite)),
                mad_log2fc=float(np.median(np.abs(finite-np.median(finite)))),
                p10_log2fc=float(np.quantile(finite,.1)),p90_log2fc=float(np.quantile(finite,.9))))
            sampling.append(dict(guide=row.guide,repeat=repeat+1,
                                 source_rows=obs.iloc[chosen].source_row.astype(int).tolist()))
            if repeat==0:matched_first.append(effect)
    if not records:raise ValueError('No guide reaches declared descriptive coverage')
    return pd.DataFrame(records),np.array(means),np.array(reference_means),np.array(effects),np.array(matched_first),pd.DataFrame(metric_rows),sampling


def pair_effects(a,b):
    keep=np.isfinite(a)&np.isfinite(b)
    x=np.asarray(a)[keep];y=np.asarray(b)[keep]
    valid=len(x)>1 and np.std(x)>1e-12 and np.std(y)>1e-12
    active=(np.abs(x)>1e-12)&(np.abs(y)>1e-12)
    return dict(features=int(keep.sum()),pearson=float(np.corrcoef(x,y)[0,1]) if valid else None,
                mean_absolute_difference=float(np.mean(np.abs(x-y))) if len(x) else None,
                sign_agreement_nonzero=float(np.mean(np.sign(x[active])==np.sign(y[active]))) if active.any() else None,
                nonzero_sign_features=int(active.sum()),correlation_defined=bool(valid))


def target_support(normalized, obs, overlay, records, effects, reference, pseudocount=1.):
    """Finite per-target support, including absent and single-guide targets."""
    rows=[];pairs=[];guides=obs.guide_id.astype(str).to_numpy()
    lookup={g:i for i,g in enumerate(records.guide)}
    center=np.asarray(normalized[reference].mean(0)).ravel()
    for target,library in overlay.loc[~overlay.source_control_type.isin(TYPES[:2])].groupby('target_gene_name',sort=True):
        library=library.sort_values('guide')
        members=records.loc[records.guide.isin(library.guide)].sort_values(['cells','guide'],ascending=[False,True])
        target_ix=np.flatnonzero(np.isin(guides,library.guide))
        largest=members.iloc[0] if len(members) else None
        corr=[]
        for a,b in combinations(members.itertuples(),2):
            result=pair_effects(effects[lookup[a.guide]],effects[lookup[b.guide]])
            if result['correlation_defined']:corr.append(result['pearson'])
            pairs.append(dict(target=target,guide_a=a.guide,guide_b=b.guide,cells_a=a.cells,cells_b=b.cells,**result))
        observed=pd.Series(guides[target_ix]).value_counts()
        largest_observed=str(observed.sort_values(ascending=False,kind='stable').index[0]) if len(observed) else None
        rest=target_ix[guides[target_ix]!=largest_observed] if largest_observed is not None else np.array([],int)
        leave_out=None
        if len(rest) and len(target_ix):
            full=expression_effect(np.asarray(normalized[target_ix].mean(0)).ravel(),center,pseudocount)
            other=expression_effect(np.asarray(normalized[rest].mean(0)).ravel(),center,pseudocount)
            direct=members.direct_feature_index.dropna() if len(members) else []
            if len(direct):full[int(next(iter(direct)))]=np.nan;other[int(next(iter(direct)))]=np.nan
            leave_out=pair_effects(full,other)['pearson']
        rows.append(dict(target=str(target),source_control_type=str(library.source_control_type.iloc[0]),
            library_guides=len(library),observed_guides=len(observed),cells=len(target_ix),eligible_guides=len(members),
            pairs=len(corr),median_pair_r=float(np.median(corr)) if corr else None,
            largest_guide=largest_observed,largest_guide_fraction=float(observed.max()/len(target_ix)) if len(target_ix) else None,
            cells_without_largest=len(rest),pooled_vs_without_largest_r=leave_out,
            scatter_guide_a=str(members.iloc[0].guide) if len(members)>1 else None,
            scatter_guide_b=str(members.iloc[1].guide) if len(members)>1 else None,
            assessment='descriptive_assessable' if len(members)>1 and corr else 'insufficient_two_guide_profiles'))
    return pd.DataFrame(rows),pd.DataFrame(pairs)


def qc_flow(all_obs):
    if all_obs.index.has_duplicates:raise ValueError('Nonunique source barcode')
    rows=[]
    for status,part in all_obs.groupby(all_obs.guide_status.astype(str),sort=True):
        passing=int(part.qc_pass_stage1.astype(bool).sum())
        rows.append(dict(guide_status=status,total_cells=len(part),rna_qc_pass=passing,rna_qc_not_pass=len(part)-passing))
    table=pd.DataFrame(rows)
    if table.total_cells.sum()!=len(all_obs):raise ValueError('Cell flow does not conserve source rows')
    expected=all_obs.qc_pass_stage1.astype(bool)&all_obs.guide_status.astype(str).eq('single_guide')
    if not np.array_equal(expected,all_obs.main_analysis_stage1.astype(bool)):
        raise ValueError('Frozen cell-QC policy differs')
    return table
