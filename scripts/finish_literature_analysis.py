"""Reference HDBSCAN and STRING11.5 annotation, without energy admission gates."""
from pathlib import Path
from datetime import datetime
import argparse,csv,gzip,json,os,subprocess,time,importlib.metadata as md
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from analyze_expression_structure import sha256,atomic_json
from run_chip01_literature import verify_stage
from replogle2022_statistics import strong_perturbation_mask,response_feature_union,paper_hdbscan


def read_string_clusters(root):
    source=json.loads((root/'source.json').read_text())
    for entry in source['references']:
        p=root/entry['file']
        if entry['version']!='STRING11.5' or p.stat().st_size!=entry['bytes'] or sha256(p)!=entry['sha256']:
            raise ValueError('STRING reference version/hash mismatch')
    proteins={}
    with gzip.open(root/'string11.5_protein_info.gz','rt',encoding='utf8') as stream:
        reader=csv.DictReader(stream,delimiter='\t');id_key=reader.fieldnames[0]
        for row in reader:proteins[row[id_key]]=row['preferred_name']
    descriptions={}
    with gzip.open(root/'string11.5_clusters_info.gz','rt',encoding='utf8') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):descriptions[row['cluster_id']]=row['best_described_by']
    members={}
    with gzip.open(root/'string11.5_clusters_proteins.gz','rt',encoding='utf8') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            gene=proteins.get(row['protein_id'])
            if gene:members.setdefault(row['cluster_id'],set()).add(gene)
    return members,descriptions,set(proteins.values()),source


def database_support(group_targets,members,descriptions):
    records=[];group=set(group_targets)
    for key,genes in sorted(members.items()):
        overlap=group&genes
        if genes and len(overlap)/len(genes)>=.75:
            records.append(dict(database='STRING11.5',database_cluster=key,description=descriptions.get(key,''),
                matched_targets=';'.join(sorted(overlap)),matched_members=len(overlap),database_members=len(genes),
                fraction_of_database_cluster=len(overlap)/len(genes),interpretation='co-response association, not physical membership validation'))
    return records


def run(args):
    if args.output.exists():raise FileExistsError('New functional output required')
    effect_state=verify_stage(args.effects);core_state=verify_stage(args.core)
    config=json.loads(args.config.read_text())
    if config['energy_significance_required_for_clustering']:raise ValueError('Forbidden additional energy gate')
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**extra):atomic_json(args.output/'status.json',dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),code_commit=commit,**extra))
    try:
        state('verify_and_cluster_reference_profiles')
        summary=pd.read_csv(args.effects/'de_summary.csv')
        energy=pd.read_csv(args.effects/'energy_results.csv')
        effects_summary=json.loads((args.effects/'summary.json').read_text())
        if not effects_summary['scientific_result']:raise ValueError('Engineering smoke cannot supply functional results')
        profiles=np.load(args.effects/'response_profiles.npz')
        if list(profiles['targets'])!=summary.target.tolist():raise ValueError('Profile target identity mismatch')
        genes=pd.read_csv(args.effects/'gene_statistics.csv.gz',index_col=0)
        if len(genes)!=profiles['values'].shape[1]:raise ValueError('Profile gene identity mismatch')
        eligible=strong_perturbation_mask(summary.cells,summary.ad_de_genes)
        features=response_feature_union(profiles['ad_statistics'],eligible,profiles['raw_mean'],profiles['raw_variance'],profiles['z_defined'])
        with threadpool_limits(config['threads']):labels,order,fit_status=paper_hdbscan(profiles['values'],eligible,features)
        summary['paper_strong_perturbation']=eligible
        summary['response_group']=['unassigned' if label<0 else 'E'+str(label+1).zfill(2) for label in labels]
        summary['energy_permutation_p']=summary.target.map(energy.set_index('target').permutation_p)
        summary['energy_bh_q']=summary.target.map(energy.set_index('target').bh_q)
        summary.to_csv(args.output/'perturbation_groups.csv',index=False)
        genes.iloc[features].to_csv(args.output/'clustering_features.csv.gz')
        state('reference_database_annotation',strong_targets=int(eligible.sum()),fit_status=fit_status)
        members,descriptions,known,annotation_source=read_string_clusters(args.annotation)
        support=[];group_rows=[];response_genes=[]
        names=genes.gene_name.astype('string').fillna('').to_numpy(dtype=str)
        for label in sorted(set(labels[labels>=0])):
            group='E'+str(label+1).zfill(2);picked=np.flatnonzero(labels==label);targets=summary.target.to_numpy()[picked]
            terms=database_support(targets,members,descriptions)
            for term in terms:term['response_group']=group;support.append(term)
            mean=profiles['values'][picked].mean(0)
            columns=features[np.argsort(np.abs(mean[features]),kind='stable')[-min(20,len(features)):]][::-1]
            for column in columns:response_genes.append(dict(response_group=group,gene_id=genes.index[column],gene_name=names[column],mean_response_z=float(mean[column])))
            group_rows.append(dict(response_group=group,targets=';'.join(targets),members=len(targets),
                known_string_target_symbols=sum(t in known for t in targets),supported_string_clusters=len(terms),
                interpretation='database-supported response association' if terms else 'unannotated response group, no invented functional label'))
        pd.DataFrame(support,columns=['response_group','database','database_cluster','description','matched_targets','matched_members','database_members','fraction_of_database_cluster','interpretation']).to_csv(args.output/'string_cluster_support.csv',index=False)
        pd.DataFrame(group_rows,columns=['response_group','targets','members','known_string_target_symbols','supported_string_clusters','interpretation']).to_csv(args.output/'group_summary.csv',index=False)
        pd.DataFrame(response_genes,columns=['response_group','gene_id','gene_name','mean_response_z']).to_csv(args.output/'group_response_genes.csv',index=False)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        core=pd.read_csv(args.core/'core_control_guides.csv')
        n_core_guides=int(core.core_control.sum())
        fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
        finite=core.mean_ad_de_genes.notna()
        ax.scatter(core.loc[finite,'cells'],core.loc[finite,'mean_ad_de_genes'],c=np.where(core.loc[finite,'core_control'],'tab:blue','tab:orange'))
        ax.axhline(config['core_control_mean_ad_de_genes_lt'],color='black',linestyle='--')
        ax.set(xlabel='Cells per NTC guide',ylabel='Mean pairwise AD-BH DE genes',title='Paper core-control selection, not a cell-quality boundary')
        fig.savefig(args.output/'core_control_selection.png',dpi=160);plt.close(fig)
        if len(order)>1 and len(features)>1:
            correlation=np.corrcoef(profiles['values'][order][:,features])
            fig,ax=plt.subplots(figsize=(10,9),layout='constrained')
            image=ax.imshow(correlation,vmin=-1,vmax=1,cmap='RdBu_r')
            ax.set(title='Strong-perturbation mean-z response correlation; HDBSCAN labels in table',xlabel='Ordered perturbations',ylabel='Ordered perturbations')
            fig.colorbar(image,ax=ax);fig.savefig(args.output/'strong_response_correlation.png',dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
        axes[0].scatter(energy.cells,energy.energy,s=12,alpha=.6);axes[0].set(xlabel='Cells per perturbation',ylabel='Full-control energy V statistic')
        axes[1].scatter(energy.permutation_p,energy.bh_q,s=12,alpha=.6);axes[1].axhline(.05,color='black',linestyle='--')
        axes[1].set(xlabel='Raw permutation p',ylabel='Across-target BH q',title='Energy values do not gate HDBSCAN')
        fig.savefig(args.output/'energy_overview.png',dpi=160);plt.close(fig)
        n_groups=len(group_rows)
        methods=dict(primary_reference=config['reference'],config=config,effects_commit=effect_state['code_commit'],core_commit=core_state['code_commit'],
            functional_commit=commit,hdbscan_version=md.version('hdbscan'),annotation=annotation_source,
            clustering_fit_status=fit_status,energy_not_a_clustering_gate=True,
            crispri_rna_30percent_filter_applied=False,dna_editing_efficiency='not_assessed',independent_biological_replicates=False,
            custom_orf_host_mapping='not_inferred',corum3_annotation_available=annotation_source['corum3_available'],
            gene_expression_program_pyMDE_branch='not part of this perturbation-cluster deliverable',exact_cross_assay_replication=False)
        atomic_json(args.output/'methods_and_scope.json',methods)
        result=dict(analysis_cells=effects_summary['analysis_cells'],core_control_cells=effects_summary['core_control_cells'],
            ad_targets=len(summary),mw_targets=int(summary.mw_tested.sum()),tested_rna_features=len(genes),
            strong_targets=int(eligible.sum()),clustering_features=len(features),hdbscan_groups=n_groups,
            strong_targets_unassigned=int((eligible&(labels<0)).sum()),fit_status=fit_status,
            energy_raw_p_lt_005=int((energy.permutation_p<.05).sum()),energy_bh_q_lt_005=int((energy.bh_q<.05).sum()),
            string_supported_groups=sum(r['supported_string_clusters']>0 for r in group_rows),
            annotation_limitation='CORUM3.0 original download unavailable; STRING11.5 exact version used; no arbitrary newer CORUM substitution')
        atomic_json(args.output/'summary.json',result)
        lines=['# 第一张芯片：文献对齐重做报告（单芯片探索性）','',
            '## 本次实际完成的内容','',
            f"沿用原始计数和原RNA QC/单sgRNA政策；按2022正文core NTC规则选{n_core_guides}条guide、{result['core_control_cells']}对照细胞。本次主分析{result['analysis_cells']}细胞，源记录未删除。",
            f"AD覆盖{result['ad_targets']}个扰动、MW覆盖{result['mw_targets']}个至少10细胞扰动，共{result['tested_rna_features']}个RNA特征。AD p值使用R kSamples尾部校准，不用SciPy .001下限；MW未测试者明确标记，不把占位p=1当实际检验。",
            f"完整核心对照的10000次标签置换能量检验：原始p<.05为{result['energy_raw_p_lt_005']}，跨扰动BH q<.05为{result['energy_bh_q_lt_005']}。两种口径分开，不据此过滤主功能聚类。",
            f"按至少25细胞且至少50 AD-BH差异基因，入选{result['strong_targets']}个强扰动；{result['clustering_features']}个特征，HDBSCAN状态{fit_status}，形成{n_groups}组，另有{result['strong_targets_unassigned']}个强扰动未成组。",
            f"其中{result['string_supported_groups']}组有STRING11.5数据库成员覆盖>=75%的注释线索。这是表达响应关联，不证明ORF编码产物属于蛋白复合物。",'',
            '## 与旧探索版本的关键区别','',
            '- core NTC采用原文成对AD而非诊断PCA异常标签；全部源细胞/guide仍保留。',
            '- 能量使用全部可用核心对照（最多5000），不是256上限；完整RNA基因AD+MW，不加绝对log2FC展示阈值作为AD差异基因定义。',
            '- 强扰动筛选不要求energy显著，HDBSCAN不被层次聚类切线替代；直接靶基因只在无歧义映射时置零，未知ORF不猜宿主。',
            '- 若入选少于4或相关谱不可定义，标记未能评价主分组；若HDBSCAN均噪声，如实报告，不改门槛求阳性。','',
            '## 仍然存在的适配和限制','',
            '- 本实验为CRISPRko，不使用CRISPRi的RNA下降30%来证明编辑成功；DNA indel/frameshift仍未评估。',
            '- 复用本批QC阈值/可用chip批次元数据，不能称所有数值与K562原文相同；方差筛选采用原始UMI方差，是对原文未明确归一化限定的一项显式解释。',
            '- 原CORUM3.0端点不可用，本轮没有拿新版数据库替代；STRING11.5原版本来源/hash见methods_and_scope.json。未注释的组保留编号，不虚构名称。',
            '- 没有独立生物学重复，不把细胞当实验重复。2026发表全文尚未完整核实，本轮不称复现该文；未执行论文另一条基因表达程序/PyMDE支线。','',
            '## 图表和完整结果','',
            '- [核心对照筛选图](core_control_selection.png)',
            '- [整体效应概览](energy_overview.png)',
            '- [扰动分组表](perturbation_groups.csv)',
            '- [分组数据库支持](string_cluster_support.csv)',
            '- [分组平均响应的高幅度基因](group_response_genes.csv)',
            '- [完整DE索引](../effects/de_summary.csv)',
            '- [AD候选基因摘要](../effects/top_ad_genes.csv.gz)',
            '- [完整能量检验表](../effects/energy_results.csv)']
        if (args.output/'strong_response_correlation.png').exists():lines.append('- [强扰动响应相关图](strong_response_correlation.png)')
        (args.output/'07_literature_redo_report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
        atomic_json(args.output/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in args.output.iterdir() if p.is_file() and p.name not in ['status.json','owner.lock','output_manifest.json']})
        state('complete',elapsed_seconds=time.monotonic()-start,summary=result,source_modified=False)
    except Exception as error:state('failed',error=repr(error));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['effects','core','annotation','config','output']:p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
