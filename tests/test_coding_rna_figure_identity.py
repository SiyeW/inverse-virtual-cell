from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from plot_coding_rna_states import save,plt


def test_dotted_transcript_ids_do_not_collapse_or_break_gallery_paths(tmp_path):
    identifiers=['PB.26158.49_chr2_1_2485_32_236_uORF_ATG',
                 'PB.26158.50_chr2_2_2500_90_270_uORF_ATG']
    for identifier in identifiers:
        fig,ax=plt.subplots(figsize=(2,1));ax.plot([0,1],[0,1])
        save(fig,tmp_path/identifier,svg=True)
    assert {p.name for p in tmp_path.iterdir()}=={identifier+'.'+suffix for identifier in identifiers for suffix in ('png','svg')}
