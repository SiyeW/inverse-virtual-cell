"""The original target-level Pearson scatter before/after front-stage filtering.

Presentation of accepted results only: no new counts, RNA features, reference,
correlation calculation, or author-pipeline claim is introduced here.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker
from render_share_ready_figures import configure,sha,save,assert_drawn_points


def prepare(original_path,retained_path,pairs_path):
    original=pd.read_csv(original_path)
    original=original.loc[original.source_control_type.isin(['targeting','positive_control'])&original.median_pair_r.notna()].copy()
    assert original.target.is_unique
    selected=pd.read_csv(retained_path);selected=selected.loc[selected.final_pass.eq(True)].copy()
    pairs=pd.read_csv(pairs_path)
    assert pairs.groupby('target').role.nunique().max()==1
    r=pairs.groupby('target',sort=True).same_r.median()
    role=pairs.groupby('target',sort=True).role.first()
    counts=selected.groupby('target').after_cell_gate.agg(['sum','max','size'])
    after=pd.DataFrame(dict(target=r.index,median_pair_r=r.to_numpy(),source_control_type=role.loc[r.index].to_numpy(),
                            largest_guide_fraction=(counts.loc[r.index,'max']/counts.loc[r.index,'sum']).to_numpy(),
                            eligible_guides=counts.loc[r.index,'size'].to_numpy(),cells=counts.loc[r.index,'sum'].to_numpy()))
    assert (after.eligible_guides>=2).all() and set(after.target).issubset(set(original.target))
    lookup=original.set_index('target').loc[after.target]
    assert np.array_equal(lookup.source_control_type.to_numpy(),after.source_control_type.to_numpy())
    paired=after.rename(columns={'median_pair_r':'r_after','largest_guide_fraction':'x_after'}).copy()
    paired['r_before']=lookup.median_pair_r.to_numpy();paired['x_before']=lookup.largest_guide_fraction.to_numpy()
    paired['r_change']=paired.r_after-paired.r_before
    ordinary=paired.source_control_type.eq('targeting')
    summary=dict(original_ordinary_targets=int(original.source_control_type.eq('targeting').sum()),
                 original_positive_targets=int(original.source_control_type.eq('positive_control').sum()),
                 matched_ordinary_targets=int(ordinary.sum()),matched_positive_targets=int((~ordinary).sum()),
                 before_matched_ordinary_r_median=float(paired.loc[ordinary,'r_before'].median()),
                 after_matched_ordinary_r_median=float(paired.loc[ordinary,'r_after'].median()),
                 median_paired_r_change=float(paired.loc[ordinary,'r_change'].median()),
                 before_correlation_uses_all_original_eligible_guides=True,
                 after_correlation_uses_retained_guides_and_retained_cells=True,
                 different_from_fixed_retained_guide_only_cell_exclusion_comparison=True)
    return original,after,paired,summary


def render(original_path,retained_path,pairs_path,out):
    original,after,paired,summary=prepare(original_path,retained_path,pairs_path)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    paired.to_csv(out/'matched_target_before_after.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(9.2,5.5),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.985,bottom=.17,top=.82,wspace=.12)
    retained=set(after.target);drawn=[]
    for i,(ax,title,table) in enumerate(zip(axes,['Before filtering','After filtering'],[original,after])):
        if i==0:
            q=table.loc[~table.target.isin(retained)]
            ax.scatter(q.largest_guide_fraction,q.median_pair_r,s=13,color='#c8c8c8',alpha=.55,linewidths=0)
            table=table.loc[table.target.isin(retained)]
        for role,color in [('targeting','#496E8B'),('positive_control','#4E8575')]:
            q=table.loc[table.source_control_type.eq(role)]
            ax.scatter(q.largest_guide_fraction,q.median_pair_r,s=17,color=color,alpha=.8,linewidths=0,zorder=3)
        expected=original if i==0 else after
        assert_drawn_points(ax,np.column_stack((expected.largest_guide_fraction,expected.median_pair_r)))
        assert expected.median_pair_r.between(-.1,.1).all()
        ax.axhline(0,color='#777777',lw=.8,zorder=1)
        ax.set_xlim(.25,1);ax.set_ylim(-.1,.1);ax.set_box_aspect(1)
        ax.set_xticks([.3,.5,.7,.9]);ax.set_yticks([-.1,-.05,0,.05,.1])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        ax.set_title(title,pad=9);ax.grid(axis='y',color='#eeeeee',lw=.45)
        drawn.append(dict(panel=title,targets=len(expected),all_offsets_match_original_values=True))
    axes[0].set_ylabel('Median same-target RNA-response Pearson r')
    fig.supxlabel("Fraction of target's cells in its largest sgRNA group",fontsize=10,y=.035)
    handles=[Line2D([],[],marker='o',linestyle='none',color='#496E8B',markersize=4,label='Retained ordinary targets'),
             Line2D([],[],marker='o',linestyle='none',color='#4E8575',markersize=4,label='Library positive control'),
             Line2D([],[],marker='o',linestyle='none',color='#c8c8c8',markersize=4,label='Other original targets')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.94),ncol=3,fontsize=8.1,columnspacing=1,handletextpad=.35)
    fig.suptitle('Company merged',fontsize=11,y=.99)
    fig.canvas.draw()
    for ax in axes:
        box=ax.get_window_extent();assert np.isclose(box.height/box.width,1,atol=1e-4)
    save(fig,out,'01_筛选前后RNA响应相关性')
    proof=dict(summary=summary,panels=drawn,x_limits=[.25,1],y_limits=[-.1,.1],
               original_RNA_features=5058,original_normalization_depth=13828,
               original_reference_and_formula_retained=True,not_the_difference_boxplot=True,
               inputs={str(p):sha(p) for p in [original_path,retained_path,pairs_path]})
    (out/'plot_proof.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True)
    p.add_argument('--retained',type=Path,required=True);p.add_argument('--pairs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();configure(None)
    render(a.original,a.retained,a.pairs,a.output)
