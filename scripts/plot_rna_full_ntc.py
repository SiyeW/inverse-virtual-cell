"""Single publication-scale RNA-state exploration panel, no guide selection."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
from explore_rna_full_ntc import energy_u


def draw(folder,output,identity_table,guide_design=None,highlight_config=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter,MaxNLocator
    folder=Path(folder);out=Path(output);out.mkdir(exist_ok=True)
    f=pd.read_csv(folder/'RNA_state_all_records.csv');s=json.loads((folder/'summary.json').read_text())
    old=pd.read_csv(identity_table).set_index('guide');assert len(f)==785 and f.guide.is_unique and int(f.query_records.sum())==71664
    assert set(f.guide)==set(old.index)
    for r in f.itertuples():
        assert r.query_records==int(old.loc[r.guide,'n_each'])
        assert r.reference_records==19342-(r.query_records if r.role=='non_targeting_control' else 0)
        assert np.isclose(energy_u(r.query_records,r.reference_records,r.within_query_sum,r.cross_sum,r.within_reference_sum),r.energy_unbiased,atol=1e-10)
    assert not s['query_cells_downsampled'] and not s['reference_NTC_downsampled']
    target=f.role.eq('targeting');positive=f.role.eq('positive_control');ntc=f.role.eq('non_targeting_control')
    assert target.sum()==659 and positive.sum()==19 and ntc.sum()==106
    assert not (target & positive).any() and int((target | positive | ntc).sum())==784
    highlighted=f.guide.isin([]);highlight_summary=[];highlight_points=[];not_plotted=[];offsets={}
    assert (guide_design is None)==(highlight_config is None)
    if highlight_config is not None:
        design=pd.read_csv(guide_design,keep_default_na=False);cfg=json.loads(Path(highlight_config).read_text())
        assert design.guide.is_unique and set(f.guide)<=set(design.guide)
        genes={x['gene_id']:x['gene_name'] for x in cfg['genes']}
        assert len(genes)==len(cfg['genes']) and set(genes)<=set(design.target_gene_id)
        chosen=design[design.target_gene_id.isin(genes)].copy()
        assert chosen.source_control_type.eq('targeting').all()
        highlighted=f.guide.isin(chosen.guide)
        assert (highlighted & target).equals(highlighted)
        offsets=cfg.get('label_offsets_points',{})
        for gene_id,gene_name in genes.items():
            rows=chosen[chosen.target_gene_id.eq(gene_id)];shown=rows.guide.isin(f.guide)
            highlight_summary.append(dict(gene_id=gene_id,gene_name=gene_name,designed_guides=len(rows),plotted_guides=int(shown.sum())))
        for row in chosen.itertuples():
            if row.guide not in set(f.guide):
                n=int(row.single_guide_after_stage1_qc)
                assert n<10
                not_plotted.append(dict(guide=row.guide,gene_id=row.target_gene_id,gene_name=genes[row.target_gene_id],records=n,
                    reason='no_QC_single_guide_records' if n==0 else 'below_original_10_record_minimum'))
        for row in f[highlighted].itertuples():
            entry=chosen.set_index('guide').loc[row.guide]
            highlight_points.append(dict(guide=row.guide,gene_id=entry.target_gene_id,gene_name=genes[entry.target_gene_id],
                records=int(row.query_records),energy_unbiased=float(row.energy_unbiased)))
        assert len(highlight_points)+len(not_plotted)==len(chosen)
    other_target=target & ~highlighted
    plt.rcParams.update({'font.sans-serif':['Arial','DejaVu Sans'],'font.size':10,'axes.unicode_minus':False,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(6.7,4.9));fig.subplots_adjust(left=.14,right=.97,top=.90,bottom=.17)
    ax.axhline(0,color='#8a8d91',lw=.65,zorder=1)
    ax.scatter(f.loc[other_target,'query_records'],f.loc[other_target,'energy_unbiased'],s=12,c='#bd8050',alpha=.48,linewidths=0,label=f'Other targeting sgRNAs (n = {int(other_target.sum())})',zorder=2)
    ax.scatter(f.loc[ntc,'query_records'],f.loc[ntc,'energy_unbiased'],s=16,c='#377394',alpha=.82,linewidths=.25,edgecolors='white',label='NTC sgRNAs (n = 106)',zorder=3)
    ax.scatter(f.loc[positive,'query_records'],f.loc[positive,'energy_unbiased'],s=25,c='#e60012',alpha=.95,linewidths=.35,edgecolors='white',label='Positive controls in library design (n = 19)',zorder=4)
    if highlighted.any():
        ax.scatter(f.loc[highlighted,'query_records'],f.loc[highlighted,'energy_unbiased'],s=30,marker='D',c='#703c99',alpha=1,linewidths=.4,edgecolors='white',label=f'Highlighted sgRNAs (n = {int(highlighted.sum())})',zorder=5)
        for point in highlight_points:
            dx,dy=offsets.get(point['guide'],[6,9])
            ax.annotate(point['guide'],xy=(point['records'],point['energy_unbiased']),xytext=(dx,dy),textcoords='offset points',
                ha='left' if dx>=0 else 'right',va='bottom' if dy>=0 else 'top',fontsize=7,color='#4d286b',zorder=6,
                arrowprops=dict(arrowstyle='-',color='#703c99',lw=.4,alpha=.7,shrinkA=1,shrinkB=3))
    ax.set_xscale('log');ticks=[10,20,50,100,200,500];ax.set_xticks(ticks);ax.xaxis.set_major_formatter(ScalarFormatter());ax.minorticks_off()
    ax.set_xlim(9.3,max(f.query_records.max(),500)*1.12)
    ax.set_xlabel('Cell barcodes per sgRNA',labelpad=8)
    ax.set_ylabel('RNA-state difference from NTC (E-distance estimate)',labelpad=8)
    ax.set_title('RNA-state differences (one point per sgRNA)',loc='left',fontsize=12,pad=12)
    ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_linewidth(.8);ax.tick_params(width=.7,length=3.5);ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.legend(frameon=False,fontsize=8.5,loc='upper right',handlelength=1.5,labelspacing=.6)
    for ext in ['png','svg']:fig.savefig(out/f'整体RNA状态差异_全量NTC.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    description=dict(all785_groups_checked=True,all71664_records_retained=True,full_NTC_reference=True,no_reference_sampling=True,no_response_filter_or_hit_list=True,
        plotted_targeting_groups=int(other_target.sum()),plotted_highlighted_groups=int(highlighted.sum()),plotted_positive_control_groups=int(positive.sum()),
        plotted_targeting_total_including_design_controls=int((target | positive).sum()),
        plotted_NTC_groups=int(ntc.sum()),intergenic_group_retained_in_table_not_target_vs_NTC_comparison=True,
        positive_control_identity='pre-existing library design annotation, not experimentally established successful KO',
        all_eligible_design_positive_guides_shown=True,
        positive_control_color='#e60012',highlighted_genes=highlight_summary,highlighted_points=highlight_points,
        highlighted_guides_not_plotted=not_plotted,all_eligible_requested_guides_shown=True,
        guide_design_SHA256=hashlib.sha256(Path(guide_design).read_bytes()).hexdigest() if guide_design is not None else None,
        highlight_config_SHA256=hashlib.sha256(Path(highlight_config).read_bytes()).hexdigest() if highlight_config is not None else None,
        plotted_coordinate_SHA256=hashlib.sha256(f.loc[target | positive | ntc,['guide','query_records','energy_unbiased']].sort_values('guide').to_csv(index=False).encode('utf-8')).hexdigest(),
        source_table_SHA256=hashlib.sha256((folder/'RNA_state_all_records.csv').read_bytes()).hexdigest(),
        observation_unit='one sgRNA and all eligible records assigned to that sgRNA; not one target gene or one cell',
        model_reference_band_shown=False,confidence_interval_shown=False,figure_language='English')
    (out/'plot_proof.json').write_text(json.dumps(description,indent=2));print(json.dumps(description))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--identity-table',type=Path,required=True)
    p.add_argument('--guide-design',type=Path);p.add_argument('--highlight-config',type=Path)
    a=p.parse_args();draw(a.folder,a.output,a.identity_table,a.guide_design,a.highlight_config)
