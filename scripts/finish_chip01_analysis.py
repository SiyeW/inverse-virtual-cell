"""Offline pathway ORA, perturbation-profile clustering and Chinese report.

Uses the tested expressed-gene universe, never uploads laboratory gene lists.
All associations are exploratory within one chip; clusters are not proof of
shared function and custom target names are not mapped to guessed host genes.
"""
from pathlib import Path
from datetime import datetime
import argparse,csv,gzip,json,os,re,subprocess,time
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster,leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import hypergeom,spearmanr
from threadpoolctl import threadpool_limits
from analyze_expression_structure import sha256,atomic_json
from perturbation_statistics import bh_adjust


def read_gmt(path):
    pathways={}
    for line in Path(path).read_text(encoding='utf8').splitlines():
        fields=line.split('\t')
        if len(fields)<3:continue
        match=re.search(r'R-HSA-\d+',fields[0]+' '+fields[1])
        if not match:continue
        key=match.group();genes=set(g.strip() for g in fields[2:] if g.strip())
        if key in pathways:raise ValueError('Duplicate human pathway identifier')
        pathways[key]=(fields[0],genes)
    if not pathways:raise ValueError('No human R-HSA pathways parsed')
    return pathways


def overrepresentation(query,background,pathways,min_size=10,max_size=500):
    """One-sided hypergeometric ORA; BH across ALL eligible pathways.

    Background is the common tested gene-symbol universe, not the genome.
    Gene IDs mapping to the same symbol collapse once; query/background are sets.
    """
    background=set(background);query=set(query)&background
    if not background:raise ValueError('Empty gene universe')
    records=[]
    for key,(name,genes) in sorted(pathways.items()):
        genes=genes&background
        if not min_size<=len(genes)<=max_size:continue
        hits=genes&query
        records.append(dict(pathway_id=key,pathway_name=name,pathway_genes_in_background=len(genes),
          query_genes=len(query),background_genes=len(background),overlap=len(hits),genes=';'.join(sorted(hits)),
          pvalue=float(hypergeom.sf(len(hits)-1,len(background),len(genes),len(query))) if query else 1.0))
    table=pd.DataFrame(records)
    if len(table):table['bh_q']=bh_adjust(table.pvalue.to_numpy())
    return table


def profile_clusters(values,eligible,max_genes=2000,distance_cut=.7,min_members=3):
    values=np.asarray(values,dtype=float);eligible=np.asarray(eligible,dtype=bool)
    if values.ndim!=2 or eligible.shape!=(len(values),) or not np.isfinite(values).all():raise ValueError('Invalid profiles')
    indices=np.flatnonzero(eligible & (values.std(1)>1e-8))
    groups=np.full(len(values),-1,dtype=int)
    if len(indices)<min_members:return groups,indices,np.array([],dtype=int)
    variance=values[indices].var(0)
    variable=np.flatnonzero(variance>1e-12)
    columns=variable[np.argsort(variance[variable],kind='stable')[-min(max_genes,len(variable)):]]
    if len(columns)<2:return groups,indices,np.array([],dtype=int)
    distances=pdist(values[indices][:,columns],metric='correlation')
    if not np.isfinite(distances).all():raise ValueError('Nonfinite profile correlation')
    tree=linkage(np.clip(distances,0,2),method='average')
    provisional=fcluster(tree,t=distance_cut,criterion='distance')
    next_id=0
    for cluster in sorted(set(provisional)):
        members=indices[provisional==cluster]
        if len(members)>=min_members:groups[members]=next_id;next_id+=1
    return groups,indices[leaves_list(tree)],columns


def verify_stage(folder):
    state=json.loads((folder/'status.json').read_text())
    if state['stage']!='complete':raise ValueError('Dependent stage not complete')
    manifest=json.loads((folder/'output_manifest.json').read_text())
    for rel,entry in manifest.items():
        p=folder/rel
        if p.stat().st_size!=entry['bytes'] or sha256(p)!=entry['sha256']:raise ValueError('Stage output identity mismatch: '+rel)
    return state,manifest


def run(args):
    if args.output.exists():raise FileExistsError('Do not overwrite prior functional run')
    effects_state,effect_manifest=verify_stage(args.effects)
    structure_state,structure_manifest=verify_stage(args.structure)
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**kw):atomic_json(args.output/'status.json',dict(stage=stage,pid=os.getpid(),time=datetime.now().astimezone().isoformat(),code_commit=commit,**kw))
    try:
        state('verify_annotation_and_mapping')
        pathways=read_gmt(args.gmt)
        annotation=json.loads(args.gmt.with_name('source.json').read_text())
        if sha256(args.gmt)!=annotation['sha256']:raise ValueError('Annotation source SHA256 mismatch')
        universe=pd.read_csv(args.effects/'tested_gene_universe.csv',index_col=0)
        symbols=universe.gene_name.fillna('').astype(str)
        background=set(symbols)-{''}
        annotated=set.union(*(genes for _,genes in pathways.values()))&background
        summary=pd.read_csv(args.effects/'de_summary.csv')
        energy=pd.read_csv(args.effects/'energy_results.csv')
        profiles=np.load(args.effects/'perturbation_expression_profiles.npz')
        if list(profiles['targets'])!=summary.target.tolist():raise ValueError('Profile/target order mismatch')
        if profiles['values'].shape[1]!=len(universe):raise ValueError('Profile/gene universe mismatch')
        summary['target_code']=['T'+str(i).zfill(3) for i in range(len(summary))]
        summary['energy_bh_q']=summary.target.map(energy.set_index('target').bh_q)
        eligible=(summary.cells>=20)&(summary.exploratory_de_genes>=5)&(summary.energy_bh_q<.05)
        labels,order,columns=profile_clusters(profiles['values'],eligible.to_numpy())
        summary['expression_effect_group']=['unassigned' if i<0 else 'E'+str(i+1).zfill(2) for i in labels]
        summary.to_csv(args.output/'perturbation_groups.csv',index=False)
        universe.iloc[columns].to_csv(args.output/'profile_clustering_features.csv')
        matrix=profiles['values']
        if len(order)>1 and len(columns)>1:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            correlation=np.corrcoef(matrix[order][:,columns])
            fig,ax=plt.subplots(figsize=(12,10),layout='constrained')
            image=ax.imshow(correlation,vmin=-1,vmax=1,cmap='RdBu_r')
            positions=np.arange(len(order)) if len(order)<=60 else np.arange(0,len(order),max(1,len(order)//40))
            names=summary.target_code.to_numpy()[order][positions]
            ax.set_xticks(positions,names,rotation=90,fontsize=6);ax.set_yticks(positions,names,fontsize=6)
            ax.set_title('Perturbation expression-effect correlation (exploratory)');fig.colorbar(image,ax=ax)
            fig.savefig(args.output/'perturbation_effect_correlation.png',dpi=170);plt.close(fig)
        hits=[];robustness=[]
        with gzip.open(args.output/'pathway_enrichment_all_tests.csv.gz','wt',encoding='utf8',newline='') as stream:
            writer=None
            for i,row in summary.iterrows():
                state('pathway_enrichment',completed_targets=int(i),total_targets=len(summary))
                de=pd.read_csv(args.effects/row.file,index_col=0)
                if not de.index.equals(universe.index):raise ValueError('DE gene order/universe mismatch')
                primary=(de.bh_q<.05)&(de.log2_fold_change.abs()>=.25)
                secondary=(de.ntc_sensitivity_bh_q<.05)&(de.ntc_sensitivity_log2_fold_change.abs()>=.25)
                pset=set(de.index[primary]);sset=set(de.index[secondary])
                robustness.append(dict(target=row.target,primary_de=len(pset),sensitivity_de=len(sset),
                  overlap=len(pset&sset),jaccard=len(pset&sset)/len(pset|sset) if pset|sset else 1.0,
                  primary_de_fraction_retained=len(pset&sset)/len(pset) if pset else np.nan,
                  log2fc_spearman=float(spearmanr(de.log2_fold_change,de.ntc_sensitivity_log2_fold_change).statistic)))
                for direction,mask in [('up',primary & (de.log2_fold_change>0)),('down',primary & (de.log2_fold_change<0))]:
                    query=set(symbols.loc[de.index[mask]])-{''}
                    ora=overrepresentation(query,annotated,pathways)
                    if ora.empty:raise ValueError('No eligible pathways in tested universe')
                    ora.insert(0,'direction',direction);ora.insert(0,'target',row.target)
                    if writer is None:writer=csv.DictWriter(stream,fieldnames=ora.columns);writer.writeheader()
                    writer.writerows(ora.to_dict('records'))
                    good=ora[(ora.bh_q<.05)&(ora.overlap>=3)].copy()
                    if len(good):hits.append(good)
        hit_table=pd.concat(hits,ignore_index=True) if hits else pd.DataFrame(columns=['target','direction','pathway_id','pathway_name','overlap','bh_q'])
        hit_table.to_csv(args.output/'pathway_enrichment_hits.csv.gz',index=False)
        pd.DataFrame(robustness).to_csv(args.output/'control_sensitivity_summary.csv',index=False)
        # Group-level labels summarize independently tested member associations,
        # not a new p-value obtained by pretending pathways are replicates.
        if len(hit_table):
            group_map=summary.set_index('target').expression_effect_group
            descriptions=hit_table.assign(expression_effect_group=hit_table.target.map(group_map))
            descriptions=descriptions[descriptions.expression_effect_group!='unassigned']
            group_terms=descriptions.groupby(['expression_effect_group','direction','pathway_id','pathway_name']).agg(
                member_targets=('target','nunique'),minimum_member_q=('bh_q','min')).reset_index()
            group_terms.to_csv(args.output/'group_pathway_support.csv',index=False)
        atomic_json(args.output/'methods_and_scope.json',dict(
          code_commit=commit,effects_code_commit=effects_state['code_commit'],structure_code_commit=structure_state['code_commit'],
          annotation=annotation,annotation_human_pathways=len(pathways),tested_gene_ids=len(universe),background_gene_symbols=len(background),
          reactome_annotated_background_symbols=len(annotated),ora_background_gene_symbols=len(annotated),annotation_background_coverage=len(annotated)/len(background),
          pathway_test='background = tested expressed symbols intersect Reactome annotated; one-sided hypergeometric; BH over all size10–500 pathways per target and direction; display overlap>=3/q<.05',
          de_gene_selection='primary per-target gene BH q<.05 and abs(CP10k mean log2FC)>=.25; up/down separately',
          cluster='targets>=20 cells and>=5 DE genes and energy q<.05; top2000 between-target variable effects; correlation distance, average linkage, cut0.7, groups>=3',
          cluster_interpretation='descriptive expression-response similarity, not independently validated shared target function',
          biological_replication=False,dna_editing_efficiency='not_assessed',custom_target_host_mapping='not_inferred',
          outputs_are_exploratory=True,laboratory_genes_uploaded=False))
        policy=json.loads((args.structure/'qc_policy.json').read_text())
        n_groups=len(set(labels[labels>=0]))
        top=energy.sort_values('energy_full',ascending=False).head(15)
        lines=['# 第一张芯片常规分析报告（探索性）','',
          '## 先看结论','',
          f'本轮使用 {policy["main_cells"]:,} 个 RNA QC 通过且仅分配一条 sgRNA 的细胞，其中 {policy["controls"]:,} 个阴性对照；原始计数和交付不修改。',
          f'差异表达分析覆盖 {len(summary)} 个至少10细胞的扰动；能量距离覆盖 {len(energy)} 个至少20细胞的扰动。',
          f'有 {int((energy.bh_q<.05).sum())} 个扰动在本轮抽样置换检验中 q<0.05；这是单芯片内细胞级关联，不是独立生物学重复验证或 DNA 编辑成功率。',
          f'表达效应相似性得到 {n_groups} 个至少3成员的描述性分组；未成组者保留为 unassigned，不强制归类。',
          f'Reactome 解释覆盖 {len(annotated)}/{len(background)} 个测试背景基因符号；{hit_table.target.nunique() if len(hit_table) else 0} 个扰动有达到本轮阈值的通路条目。通路名称不是靶点功能被实验验证的证明。','',
          '## 图表入口','',
          '- 标准细胞结构图：../structure_v2/expression_umap.png；细胞分群不是已确认细胞类型。',
          '- 扰动效应相关图：perturbation_effect_correlation.png（若有足够信号生成）。T编号与完整靶点名称对应见 perturbation_groups.csv。',
          '- 每扰动基因变化：../effects/target_XXX_de.csv.gz；文件与靶点对应见 ../effects/de_summary.csv。',
          '- 能量距离：../effects/energy_results.csv；pathway_enrichment_hits.csv.gz 为通路候选，all_tests 保留完整检验。','',
          '## 整体表达变化靠前的扰动','',
          '| 靶点 | 细胞数 | 全细胞能量距离 | 抽样置换BH q |','|---|---:|---:|---:|']
        lines += [f'| {r.target} | {r.cells} | {r.energy_full:.4g} | {r.bh_q:.4g} |' for r in top.itertuples()]
        lines += ['', '## 解释边界与稳健性','',
          '保留全部阴性对照为主分析；排除此前8条诊断异常对照后的结果单列敏感性，未宣称它们失效。control_sensitivity_summary.csv给出基因集合重叠和效应相关；不能用对照敏感性替代残余双细胞敏感性。',
          '已有双细胞自动线缺乏可靠过滤依据，因此未硬切或宣称已清理所有双细胞。多sgRNA排除也不是完整双细胞检测；主要结论仍受残余混合记录影响。',
          '低覆盖扰动的效应不稳定；未覆盖/不显著不等于无功能。自定义ORF靶点不猜宿主基因，功能解释来自标准人类基因的响应。',
          '能量表示和抽样方案是对2022方法的明确适配，不称严格复现；未直接核对2026全文，不称其复现。细胞数不是生物学重复数。',
          '方法、参数、版本及注释来源见 methods_and_scope.json、../effects/statistics_parameters.json 和 ../structure_v2/representation_parameters.json。',
          '参考：[Replogle 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/)、[Scanpy常规流程](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html)、[Enrichr注释提供方](https://maayanlab.cloud/Enrichr/)。']
        (args.output/'04_standard_analysis_report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in args.output.iterdir() if p.is_file() and p.name not in ['status.json','owner.lock','output_manifest.json']}
        atomic_json(args.output/'output_manifest.json',manifest)
        state('complete',elapsed_seconds=time.monotonic()-start,targets=len(summary),energy_targets=len(energy),
              expression_effect_groups=n_groups,report_needs_scientific_inspection=True,source_modified=False)
    except Exception as exc:state('failed',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--effects',type=Path,required=True);p.add_argument('--structure',type=Path,required=True)
    p.add_argument('--gmt',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    with threadpool_limits(2):run(p.parse_args())
