from pathlib import Path
import sys
from types import SimpleNamespace
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_guide_sequence_coverage import parse_mappings,read_observation,mapping_command,locus_coverage


def test_splice_skips_are_not_editable_sequence_coverage():
    ordinary=SimpleNamespace(reference_start=90,cigartuples=[(0,40)])
    deletion=SimpleNamespace(reference_start=90,cigartuples=[(0,18),(2,4),(0,18)])
    spliced=SimpleNamespace(reference_start=90,cigartuples=[(0,10),(3,20),(0,20)])
    assert read_observation(ordinary,100,120)==(True,True)
    assert read_observation(deletion,100,120)==(True,True)
    assert read_observation(spliced,100,120)==(False,False)


def test_exact_unique_placement_never_selects_from_multi_hits(tmp_path):
    rows=[dict(guide=n,sequence='A'*20,target_gene_name='x',source_control_type='targeting') for n in ['one','multi','absent']]
    p=tmp_path/'mapped.sam'
    p.write_text('g0000\t0\tchr1\t101\t255\t20M\t*\t0\t0\t'+ 'A'*20+'\t*\tNH:i:1\tnM:i:0\n'+
                 'g0001\t0\tchr1\t201\t3\t20M\t*\t0\t0\t'+ 'A'*20+'\t*\tNH:i:2\tnM:i:0\n'+
                 'g0002\t4\t*\t0\t0\t*\t*\t0\t0\t'+ 'A'*20+'\t*\n')
    r=parse_mappings(p,rows)
    assert r[0]['start0']==100 and r[0]['end0']==120
    assert r[1]['mapping_status']=='ambiguous_exact_placement' and r[1]['start0'] is None
    assert r[2]['mapping_status']=='not_exactly_located'


def test_coverage_counts_families_without_calling_an_edit():
    def read(cb,ub,quality=60):
        tags={'CB':cb}
        if ub is not None:tags['UB']=ub
        return SimpleNamespace(reference_start=90,cigartuples=[(0,40)],is_unmapped=False,
                               is_secondary=False,is_supplementary=False,mapping_quality=quality,
                               has_tag=lambda k:k in tags,get_tag=lambda k:tags[k])
    reads=[read('own','u1'),read('own','u1'),read('control','u2'),read('own','u3',quality=3)]
    bam=SimpleNamespace(fetch=lambda *_:iter(reads))
    row=dict(guide='A-1',chromosome='chr1',start0=100,end0=120)
    r=locus_coverage(bam,row,{'own':'A-1','control':'NT-1'},{'NT-1'},1000)
    assert r['own_guide_distinct_CB_UB_families_spanning']==1
    assert r['NTC_distinct_CB_UB_families_spanning']==1
    assert r['primary_reads_spanning_spacer']==3
    reads[:]=[read('own',None)]
    r=locus_coverage(bam,row,{'own':'A-1'},{'NT-1'},1000)
    assert r['coverage_status']=='own_assigned_reads_span_but_UMI_identity_missing'
    assert 'editing_efficiency' not in r


def test_star_short_spacer_integration(tmp_path):
    import os,subprocess,random
    star=os.environ.get('TEST_STAR_BINARY')
    if not star:pytest.skip('Native STAR integration')
    rng=random.Random(124)
    dna=''.join(rng.choices('ACGT',k=6000))
    one=dna[1000:1020];multi=dna[2000:2020]
    dna=dna[:4000]+multi+dna[4020:]
    genome=tmp_path/'genome.fa';genome.write_text('>chrToy\n'+dna+'\n')
    index=tmp_path/'index';index.mkdir()
    r=subprocess.run([star,'--runMode','genomeGenerate','--runThreadN','2','--genomeDir',str(index),
                      '--genomeFastaFiles',str(genome),'--genomeSAindexNbases','4','--genomeChrBinNbits','8',
                      '--outFileNamePrefix',str(tmp_path/'build_')],capture_output=True,text=True,timeout=40)
    assert r.returncode==0,r.stdout+r.stderr
    changed=one[:-1]+next(b for b in 'ACGT' if b!=one[-1])
    reverse=one.translate(str.maketrans('ACGT','TGCA'))[::-1]
    rows=[dict(guide=n,sequence=s,target_gene_name='x',source_control_type='targeting') for n,s in [('one',one),('multi',multi),('absent','T'*20),('one_mismatch',changed),('reverse',reverse)]]
    fasta=tmp_path/'guides.fa';fasta.write_text(''.join(f'>g{i:04d}\n{r["sequence"]}\n' for i,r in enumerate(rows)))
    prefix=tmp_path/'map_'
    r=subprocess.run(mapping_command(star,index,fasta,prefix),capture_output=True,text=True,timeout=40)
    assert r.returncode==0,r.stdout+r.stderr
    result=parse_mappings(str(prefix)+'Aligned.out.sam',rows,lambda _chr,a,b:dna[a:b])
    assert result[0]['mapping_status']=='unique_exact_genomic_placement' and result[0]['start0']==1000
    assert result[1]['mapping_status']=='ambiguous_exact_placement'
    assert result[2]['mapping_status']=='not_exactly_located'
    assert result[3]['mapping_status']=='not_exactly_located'
    assert result[4]['mapping_status']=='unique_exact_genomic_placement'
    assert result[4]['start0']==1000 and result[4]['sequence_alignment_strand']=='-'
