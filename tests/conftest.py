"""Keep historical numerical regression tests runnable after route archival.

Only tests opt into the legacy import path. The adopted runner never does.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT/'archive/20261023_legacy_chip_analysis/scripts'))
