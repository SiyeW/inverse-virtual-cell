"""Join target-feature counts and local genomic RNA observations without nesting them.

This is an observability summary, not an editing/KO classifier. Read-cap
negatives and unplaced sites stay unknown; a positive in a capped scan is a
positive lower bound. Expression presence never uses the NTC denominator.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check_membership(expression_cells,region_cells):
    for x in [expression_cells,region_cells]:
        if not x.cell_barcode.is_unique:raise ValueError('Duplicate record identity')
    a=expression_cells.set_index('cell_barcode').sort_index();b=region_cells.set_index('cell_barcode').sort_index()
    if not a.index.equals(b.index):raise ValueError('Record populations differ')
    if not np.array_equal(a.guide_id.astype(str),b.guide_id.astype(str)):raise ValueError('Per-record guide assignments differ')
    if not np.array_equal(pd.to_numeric(a.total_counts),pd.to_numeric(b.total_counts)):raise ValueError('Full-RNA totals differ')
    return len(a)

def combine(expression,regional,placements,cells):
    if not expression.guide.is_unique or not regional.guide.is_unique or not placements.guide.is_unique:
        raise ValueError('Duplicate guide identity')
    if not expression.mapping_status.eq('exact_catalogue_id').all():raise ValueError('Unresolved target-feature mapping cannot be called undetected')
    sizes=cells.groupby('guide_id').size();expected=expression.guide.map(sizes).fillna(0)
    if not np.array_equal(expression.cells,expected):raise ValueError('Guide record counts differ')
    if not set(regional.guide)<=set(expression.guide):raise ValueError('Unexpected local-RNA guide')
    rows=[];expression=expression.copy()
    region=regional.set_index('guide');mapping=placements.set_index('guide')
    for x in expression.itertuples(index=False):
        row=dict(guide=x.guide,target_gene_id=x.target_gene_id,target=x.target,records=int(x.cells),
            target_RNA_detected=None,detected_records=None,regional_RNA_status='no_assigned_records',
            region_not_assessable_reason='',observed_RNA_families=None,max_indel_support_families=None,
            scan_capped=False,positive_lower_bound_only=False)
        if x.cells==0:
            if x.guide in region.index:raise ValueError('Local reads assigned to a guide with no records')
            rows.append(row);continue
        if not np.isfinite(x.detected_cells) or not 0<=x.detected_cells<=x.cells or x.detected_cells!=round(x.detected_cells):
            raise ValueError('Invalid target-feature detection count')
        row.update(target_RNA_detected=bool(x.detected_cells>0),detected_records=int(x.detected_cells))
        if x.guide not in region.index:
            if x.guide not in mapping.index or mapping.loc[x.guide,'mapping_status']!='ambiguous_exact_placement':
                raise ValueError('Missing local result has no documented non-unique placement')
            row.update(regional_RNA_status='not_assessable',region_not_assessable_reason='non_unique_genomic_placement')
        else:
            r=region.loc[x.guide];cap=bool(r.scan_capped)
            own=int(r.own_RNA_families_observed);alt=int(r.max_alt_families)
            if own<0 or alt<0:raise ValueError('Negative RNA molecule count')
            row.update(observed_RNA_families=own,max_indel_support_families=alt,scan_capped=cap)
            # A deletion can remove all reference bases in a target interval:
            # explicit non-conflicting allele support is still positive RNA evidence.
            if own>0 or alt>0:row.update(regional_RNA_status='observed',positive_lower_bound_only=cap)
            elif cap:row.update(regional_RNA_status='not_assessable',region_not_assessable_reason='read_limit_without_positive_evidence')
            else:row['regional_RNA_status']='not_observed'
        rows.append(row)
    return pd.DataFrame(rows)

def summarize(t):
    rows=[]
    for present in [False,True]:
        q=t[t.target_RNA_detected.eq(present)]
        for region in ['not_observed','observed','not_assessable']:
            rows.append(dict(target_RNA_detected=present,regional_RNA_status=region,guides=int(q.regional_RNA_status.eq(region).sum())))
    return rows

def run(cfg):
    for p,want in cfg['input_sha256'].items():
        if sha(p)!=want:raise ValueError('Frozen input differs: '+p)
    out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=True)
    assert not (out/'joint_summary.json').exists()
    export=json.loads(Path(cfg['region_export']).read_text());summary={};tables=[]
    placement=pd.read_csv(cfg['placements'],keep_default_na=False)
    design=pd.read_csv(cfg['design'],keep_default_na=False)
    wanted=design[design.source_control_type.isin(['targeting','positive_control'])]
    assert len(wanted)==1084 and wanted.guide.is_unique
    for cohort in export['cohorts']:
        ds=cohort['dataset_id'];paths=cfg['datasets'][ds]
        e=pd.read_csv(paths['expression']);r=pd.DataFrame(cohort['guide_rows'])
        ec=pd.read_csv(paths['expression_cells']);rc=pd.read_csv(paths['region_cells'])
        record_count=check_membership(ec,rc)
        if set(e.guide)!=set(wanted.guide) or len(e)!=1084:raise ValueError('Incomplete full-library expression input')
        indexed=wanted.set_index('guide')
        if e.target_gene_id.tolist()!=indexed.loc[e.guide,'target_gene_id'].tolist():raise ValueError('Target gene design differs')
        planned={g['guide']:g for reg in cohort['candidate']['regions'] for g in reg['guides']}
        if set(planned)!=set(r.guide):raise ValueError('Region result coverage differs from the complete accepted plan')
        for x in e[e.guide.isin(planned)].itertuples(index=False):
            if x.target!=planned[x.guide]['target']:raise ValueError('Guide target labels differ')
        t=combine(e,r,placement,rc);t.insert(0,'dataset',ds);tables.append(t)
        counts=summarize(t);observed=int(t.records.gt(0).sum());no_records=int(t.records.eq(0).sum())
        assert observed==cohort['candidate']['observed_target_guides'] and len(r)==cohort['candidate']['uniquely_placed_observed']
        assert sum(v['guides'] for v in counts)+no_records==1084
        summary[ds]=dict(design_guides=1084,record_count=record_count,observed_guides=observed,no_assigned_records=no_records,
            uniquely_placed_guides=len(r),counts=counts,
            non_unique_position=int(t.region_not_assessable_reason.eq('non_unique_genomic_placement').sum()),
            limited_without_positive=int(t.region_not_assessable_reason.eq('read_limit_without_positive_evidence').sum()),
            capped_positive_lower_bounds=int(t.positive_lower_bound_only.sum()),
            detected_RNA_but_no_regional=int((t.target_RNA_detected.eq(True)&t.regional_RNA_status.eq('not_observed')).sum()),
            no_target_feature_count_but_regional=int((t.target_RNA_detected.eq(False)&t.regional_RNA_status.eq('observed')).sum()))
    joined=pd.concat(tables,ignore_index=True);joined.to_csv(out/'joint_guide_observability.csv',index=False)
    proof=dict(summary=summary,input_sha256=cfg['input_sha256'],no_BAM_or_raw_count_read=True,
        per_barcode_population_guide_and_total_RNA_verified=True,
        expression='designated target RNA feature detection in all assigned QC single-guide records; not NTC detection or a remaining ratio',
        region='at least one BQ20 observed base inside the mapped20nt spacer, or observed allele support; MAPQ30, corrected CB/UB/GX RNA families',
        region_is_not_transcript_identity_or_complete_cutsite_coverage=True,no_DNA_or_functional_KO_efficiency=True)
    (out/'joint_summary.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2));print(json.dumps(summary))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()))
