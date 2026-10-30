"""Single publication-scale RNA-state exploration panel, no guide selection."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
from explore_rna_full_ntc import energy_u


def draw(folder,output,identity_table):
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
    target=f.role.isin(['targeting','positive_control']);ntc=f.role.eq('non_targeting_control');assert target.sum()==678 and ntc.sum()==106
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'font.size':10,'axes.unicode_minus':False,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(6.7,4.9));fig.subplots_adjust(left=.14,right=.97,top=.90,bottom=.17)
    ax.axhline(0,color='#8a8d91',lw=.65,zorder=1)
    ax.scatter(f.loc[target,'query_records'],f.loc[target,'energy_unbiased'],s=12,c='#bd8050',alpha=.58,linewidths=0,label='靶向sgRNA（678条）',zorder=2)
    ax.scatter(f.loc[ntc,'query_records'],f.loc[ntc,'energy_unbiased'],s=16,c='#377394',alpha=.82,linewidths=.25,edgecolors='white',label='NTC sgRNA（106条）',zorder=3)
    ax.set_xscale('log');ticks=[10,20,50,100,200,500];ax.set_xticks(ticks);ax.xaxis.set_major_formatter(ScalarFormatter());ax.minorticks_off()
    ax.set_xlim(9.3,max(f.query_records.max(),500)*1.12)
    ax.set_xlabel('该sgRNA对应的合格细胞记录数',labelpad=8)
    ax.set_ylabel('与NTC的RNA状态差异（E-distance估计）',labelpad=8)
    ax.set_title('整体RNA状态差异（每点一条sgRNA）',loc='left',fontsize=12,pad=12)
    ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_linewidth(.8);ax.tick_params(width=.7,length=3.5);ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.legend(frameon=False,fontsize=8.5,loc='upper right',handlelength=1.5,labelspacing=.6)
    for ext in ['png','svg']:fig.savefig(out/f'整体RNA状态差异_全量NTC.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    description=dict(all785_groups_checked=True,all71664_records_retained=True,full_NTC_reference=True,no_reference_sampling=True,no_response_filter_or_hit_list=True,
        plotted_targeting_groups=int(target.sum()),plotted_NTC_groups=int(ntc.sum()),intergenic_group_retained_in_table_not_target_vs_NTC_comparison=True,
        source_table_SHA256=hashlib.sha256((folder/'RNA_state_all_records.csv').read_bytes()).hexdigest(),
        observation_unit='one sgRNA and all eligible records assigned to that sgRNA; not one target gene or one cell',
        model_reference_band_shown=False,confidence_interval_shown=False)
    (out/'plot_proof.json').write_text(json.dumps(description,indent=2));print(json.dumps(description))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--identity-table',type=Path,required=True);a=p.parse_args();draw(a.folder,a.output,a.identity_table)
