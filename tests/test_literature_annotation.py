from pathlib import Path
import sys,gzip,json,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from finish_literature_analysis import read_string_clusters,database_support
from replogle2022_statistics import paper_hdbscan


def test_reference_database_coverage_not_group_purity_or_guessed_host():
    members={'c1':{'A','B','C','D'}};desc={'c1':'four gene complex'}
    assert database_support(['A','B','PB.A_orf'],members,desc)==[]
    result=database_support(['A','B','C','PB.unknown'],members,desc)
    assert len(result)==1 and result[0]['fraction_of_database_cluster']==.75
    assert result[0]['matched_targets']=='A;B;C'


def test_versioned_string_reader_checks_all_hashes(tmp_path):
    data={
      'string11.5_protein_info.gz':'#string_protein_id\tpreferred_name\tprotein_size\tannotation\np1\tA\t10\tx\np2\tB\t20\ty\n',
      'string11.5_clusters_info.gz':'#ncbi_taxid\tcluster_id\tcluster_size\tbest_described_by\n9606\tc1\t2\tAB complex\n',
      'string11.5_clusters_proteins.gz':'#ncbi_taxid\tcluster_id\tprotein_id\n9606\tc1\tp1\n9606\tc1\tp2\n'}
    refs=[]
    for name,text in data.items():
        p=tmp_path/name
        with gzip.open(p,'wt',encoding='utf8') as f:f.write(text)
        refs.append(dict(file=name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),version='STRING11.5'))
    (tmp_path/'source.json').write_text(json.dumps(dict(references=refs,corum3_available=False)))
    members,desc,known,source=read_string_clusters(tmp_path)
    assert members=={'c1':{'A','B'}} and desc['c1']=='AB complex' and known=={'A','B'}


def test_one_undefined_mean_profile_does_not_cancel_other_hdbscan_profiles():
    rng=np.random.default_rng(4)
    good=np.vstack([np.tile([1,1,0,0],(5,1)),np.tile([0,0,1,1],(5,1))])+rng.normal(scale=.001,size=(10,4))
    profiles=np.vstack([good,np.zeros((1,4))])
    labels,order,status=paper_hdbscan(profiles,np.ones(11,bool),np.arange(4))
    assert status=='hdbscan_evaluated_undefined_profiles_excluded' and labels[-1]==-1
    assert len(order)==10 and len(set(labels[:-1]))==2
