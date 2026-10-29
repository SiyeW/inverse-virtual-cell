from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_guided_rna_alleles import Read,raw_events,covers_reference,event_key,has_observed_base


def make(cigar,seq,start=10,flag=0):
    return Read(f'r\t{flag}\tchrT\t{start+1}\t60\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*len(seq))+'\tCB:Z:C1\tUB:Z:U1\tGX:Z:gene','G1')


def test_splice_N_does_not_become_a_deletion():
    fetch=lambda _c,a,b:'A'*(b-a)
    region=dict(guides=[dict(guide='g',start0=20,end0=40)])
    assert raw_events(make('10M20N10M','A'*20),region,fetch)==[]
    assert covers_reference(make('10M20N10M','A'*20),20,40) is None
    events=raw_events(make('10M3D10M','A'*20),region,fetch)
    assert len(events)==1 and events[0]['kind']=='DEL' and len(events[0]['ref'])-len(events[0]['alt'])==3


def test_insertion_and_reverse_SAM_sequence_semantics():
    fetch=lambda _c,a,b:'A'*(b-a)
    region=dict(guides=[dict(guide='g',start0=15,end0=35)])
    events=raw_events(make('10M2I10M','A'*10+'CG'+'A'*10,flag=16),region,fetch)
    assert len(events)==1 and events[0]['alt']=='ACG'
    assert make('20M','A'*20,flag=16).family==('G1','C1','U1','gene')
    assert covers_reference(make('20M','A'*20,flag=16),15,18)=='AAA'


def test_base_quality_and_reference_opportunity():
    fetch=lambda _c,a,b:'A'*(b-a)
    region=dict(guides=[dict(guide='g',start0=10,end0=30)])
    read=make('20M','A'*5+'C'+'A'*14)
    assert event_key(raw_events(read,region,fetch)[0])==('chrT',15,'A','C')
    read.qualities[5]=3
    assert raw_events(read,region,fetch)==[]
    assert covers_reference(read,15,16) is None
    assert covers_reference(make('10M2I10M','A'*10+'CG'+'A'*10),19,21) is None


def test_fast_opportunity_is_identical_to_per_base_definition():
    import random
    rng=random.Random(277)
    for cigar,n in [('100M',100),('10M2I20M3D30M',62),('10M100N20M',30),('3S10M1P2I20M5S',40)]:
        r=make(cigar,'A'*n)
        r.qualities=[rng.choice([3,20,40]) for _ in range(n)]
        for start in range(0,160,3):
            for length in [1,5,20]:
                expected=any(covers_reference(r,p,p+1) is not None for p in range(start,start+length))
                assert has_observed_base(r,start,start+length)==expected
