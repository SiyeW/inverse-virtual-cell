from pathlib import Path
import sys,pandas as pd,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from rna_state_overview import select_records

def cells():
    return pd.DataFrame(dict(cell_barcode=[f'cell{i}' for i in range(25)],guide_id=['A']*10+['NTC']*12+['B']*3,total_counts=[100]*25))

def test_only_guide_counts_control_selection_and_barcode_order_is_exact():
    c=cells();pc=c.cell_barcode.tolist()[::-1];selected,rows=select_records(c,pc)
    assert len(selected)==22 and selected.guide_id.value_counts().to_dict()=={'NTC':12,'A':10}
    assert [pc[i] for i in rows]==selected.cell_barcode.tolist()

def test_missing_eligible_record_is_not_silently_dropped():
    c=cells()
    with pytest.raises(ValueError,match='missing'):select_records(c,c.cell_barcode.tolist()[1:])

def test_duplicate_barcode_is_not_a_second_cell():
    c=cells();c.loc[1,'cell_barcode']=c.loc[0,'cell_barcode']
    with pytest.raises(ValueError,match='Duplicate'):select_records(c,c.cell_barcode.tolist())
