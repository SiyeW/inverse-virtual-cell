"""Put original lab/public RNA figures side by side without changing values."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

from perturbseq_core_io import atomic_json, sha256


COLORS = dict(non_targeting_control='#3c7c9f', intergenic_control='#898292',
              positive_control='#6e9d70', targeting='#cc8940')
NAMES = dict(non_targeting_control='NTC', intergenic_control='Intergenic',
             positive_control='Library-annotated positive', targeting='Targeting')


def render(paths, output, chinese=False):
    paths = list(map(Path, paths)); output = Path(output)
    data = [json.loads(p.read_text(encoding='utf8')) for p in paths]
    if len(data) < 2 or not any(d['dataset_id']=='kolf_original' for d in data):
        raise ValueError('A completed original-method KOLF result and a lab result are required')
    if len({d['dataset_id'] for d in data}) != len(data):
        raise ValueError('Duplicate comparison dataset')
    if any(d['matched_cells']!=10 or d['matched_repeats']!=8 or
           d['agreement_unit']!='one_point_per_target_median_of_guide_pairs' for d in data):
        raise ValueError('Changed method cannot stand in for the original comparison')
    output.mkdir(parents=True, exist_ok=False)
    (output/'figures').mkdir(); (output/'tables').mkdir()
    if chinese:
        font = Path('C:/Windows/Fonts/msyh.ttc')
        if not font.exists():
            raise ValueError('Requested readable Chinese font unavailable')
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    names = dict(non_targeting_control='阴性对照',intergenic_control='基因间区对照',
                 positive_control='文库标注正控',targeting='靶向sgRNA') if chinese else NAMES
    labels = dict(chip01_original='我们的首芯片',company_original='公司合并矩阵',
                  kolf_original='论文KOLF（Strong预选）') if chinese else {}
    title = lambda d: labels.get(d['dataset_id'],d['label'])
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                         'axes.unicode_minus':False, 'savefig.facecolor':'white'})
    bound = max(max(abs(np.asarray(d['density_edges']))) for d in data)
    ymax = max(max(v['density']) for d in data for v in d['roles'].values())*1.04
    rmsmax = max(max(v['rms_x']) for d in data for v in d['roles'].values())*1.04
    fig, axes = plt.subplots(2, len(data), figsize=(6*len(data),8.2), squeeze=False, layout='constrained')
    for col, d in enumerate(data):
        for role, v in d['roles'].items():
            label = f'{names[role]} ({v["guides"]})'
            axes[0,col].stairs(v['density'], d['density_edges'], color=COLORS[role], label=label)
            if v['rms_singleton']:
                axes[1,col].axvline(v['rms_x'][0], color=COLORS[role], label=label)
            else:
                axes[1,col].step(v['rms_x'], v['rms_y'], where='post', color=COLORS[role], label=label)
        for row in [0,1]:
            axes[row,col].legend(frameon=False, fontsize=8)
        axes[0,col].axvline(0, color='black', linestyle='--', linewidth=.7)
        axes[0,col].set(xlim=(-bound,bound), ylim=(0,ymax), xlabel='RNA变化（log2FC）' if chinese else 'RNA log2FC',
                        ylabel='密度' if chinese else 'Density',
                        title=title(d)+(f'\n{d["response_features"]:,}项RNA；每次10个细胞' if chinese else f'\n{d["response_features"]:,} RNA; 10 cells per draw'))
        axes[1,col].set(xlim=(0,rmsmax), ylim=(0,1.02), xlabel='RNA整体变化幅度（RMS）' if chinese else 'RMS RNA log2FC',
                        ylabel='累计比例' if chinese else 'Cumulative fraction',
                        title='每条sgRNA重复8次，取中位数' if chinese else 'Median of 8 draws for each sgRNA')
    fig.suptitle('原方法：RNA变化分布与整体变化幅度（相同显示尺度）' if chinese else
                 'Original method: RNA changes and response magnitude; shared display scales')
    fig.savefig(output/'figures/01_original_RNA_effects_comparison.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    fig, axes = plt.subplots(2, len(data), figsize=(6*len(data),8), squeeze=False, layout='constrained')
    max_count = max(max(d['target_coverage']) for d in data)*1.12
    for col, d in enumerate(data):
        axes[0,col].bar(['0','1','2','3+'], d['target_coverage'], color='#3c7c9f')
        axes[0,col].set(ylim=(0,max_count), xlabel='每个靶点的合格sgRNA数' if chinese else 'Eligible sgRNAs per target',
                        ylabel='靶点数' if chinese else 'Targets', title=title(d))
        for role, v in d['agreement'].items():
            axes[1,col].scatter(v['largest_fraction'], v['median_pair_r'], s=18, alpha=.65,
                                color=COLORS[role], label=f'{names[role]} ({len(v["targets"])})')
        axes[1,col].axhline(0, color='grey', linewidth=.8)
        axes[1,col].set(xlim=(0,1), ylim=(-1,1), xlabel='最大sgRNA占该靶点细胞的比例' if chinese else 'Largest sgRNA cell fraction',
                        ylabel='同靶点相关系数中位数' if chinese else 'Median Pearson r within target',
                        title='各sgRNA使用全部合格细胞' if chinese else 'All qualified cells per sgRNA')
        axes[1,col].legend(frameon=False, fontsize=8)
    fig.suptitle('原方法：同靶点不同sgRNA的一致性（相同显示尺度）' if chinese else
                 'Original method: agreement between sgRNAs targeting the same gene; shared scales')
    fig.savefig(output/'figures/02_original_agreement_comparison.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    rows = []
    for d in data:
        target = d['agreement'].get('targeting', {})
        rows.append(dict(dataset=title(d), dataset_id=d['dataset_id'], analysis_cells=d['analysis_cells'],
                         response_RNA=d['response_features'],
                         NTC_RMS=d['roles'].get('non_targeting_control',{}).get('rms_median'),
                         targeting_RMS=d['roles'].get('targeting',{}).get('rms_median'),
                         assessed_targeting_targets=len(target.get('targets', [])),
                         median_target_r=float(np.median(target['median_pair_r'])) if target else None))
    pd.DataFrame(rows).to_csv(output/'tables/original_comparison_facts.csv',index=False)
    atomic_json(output/'plot_provenance.json',dict(
        inputs=[dict(path=str(p), sha256=sha256(p)) for p in paths],
        scientific_values_changed=False, density_bins='original_source_specific_bins_preserved',
        display_limits_shared=True, original_RNA_rule_not_replaced_with_HVG=True,
        published_data_scope='preselected_Strong_release',
        lab_merged_scope='company_merged_records_diagnostic',
        method_validity_decided_by_renderer=False,
        interpretation_required='Reference responses, controls, precision and cohort differences; no dataset quality ranking'))
    manifest = {p.relative_to(output).as_posix():dict(bytes=p.stat().st_size,sha256=sha256(p))
                for p in output.rglob('*') if p.is_file()}
    atomic_json(output/'output_manifest.json',manifest)
    return dict(output=str(output), files=len(manifest), total_bytes=sum(v['bytes'] for v in manifest.values()))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exports',nargs='+',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--chinese',action='store_true')
    args=parser.parse_args()
    print(json.dumps(render(args.exports,args.output,args.chinese)))
