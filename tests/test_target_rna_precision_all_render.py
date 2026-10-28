from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
from render_target_rna_precision_all import whole_rows


def test_all_preexisting_blue_points_included_even_zero_and_weak_effects():
    t=pd.DataFrame(dict(guide=['z','b','a','missing','few'],target_gene_id=['g2','g1','g1','g3','g4'],cells=[50,20,20,20,3],NTC_mean_per_10000=[1,1,1,0,1],NTC_detected_fraction=[.1,.2,.2,0,.9],RNA_remaining_ratio=[0,1,4,0,1],source_control_type=['targeting','positive_control','targeting','targeting','targeting']))
    q=whole_rows(t)
    assert q.guide.tolist()==['a','b','z'] and q.plot_row.tolist()==[1,2,3]
    t['RNA_remaining_ratio']=[20,0,.001,1,1]
    assert whole_rows(t).guide.tolist()==q.guide.tolist()
