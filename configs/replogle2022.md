# Replogle 2022 aligned reanalysis contract

Primary reference: Replogle et al., Cell 2022, DOI
[10.1016/j.cell.2022.05.013](https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/),
STAR Methods: internal normalization, energy distance, AD/MW differential
expression, and global clustering of strong perturbations.

This contract is fixed before results, not tuned to make positive findings.
`replogle2022.json` records all numerical settings. It replaces the earlier
exploratory analysis's custom energy-significance clustering gate.

## Preserve versus recompute

Retain immutable input counts, barcode/guide assignments and the existing RNA
QC/single-guide policy. This is not a rerun from FASTQ and does not introduce
another doublet classifier. Recompute core-control selection, internal
normalization, energy tests, AD and MW statistics, strong-perturbation feature
selection, HDBSCAN and functional interpretation in new output directories.

## Reference decisions

- Core NTC: all candidate-guide pairs, genes with mean raw UMI >1, two-sample
  AD, BH across the genes within each pair (the author's library applies
  correction within contrasts; use the paper's BH rather than its old code
  default). Retain guides with mean DE count <8, the paper's genome-scale
  setting. Do not automatically replace this threshold after seeing results.
- UMI normalization: each cell totals the median core-control UMI; relative
  gene z scores use core controls within the available batch metadata.
- Energy: genes with mean raw UMI >0.5, upper z clip10, PCA20, at most5000
  control cells (all available if fewer), all cells of a tested perturbation,
  10000 label permutations. Report raw Monte Carlo p and across-target BH q
  separately. Do not reuse the old <=256 control cap. Memory budget failure
  is an explicit failure, never silent resampling.
- DE: test supplied RNA features with finite positive core-control SD, including
  defined-z genes not admitted to the PCA. Store zero-control-SD features but
  mark them untested with NaN p/q; they are outside the BH testing family.
  AD statistic uses SciPy-compatible midranks; tail calibration uses an
  externally pinned R kSamples `ad.pval(m=1,version=1)` grid. Off-grid reference
  relative error must be <=0.1%. MW uses two-sided asymptotic p for >=10 cells.
  BH is across all tested features per target, not per computational block.
  Informational fold changes are not an added DEG threshold.
- Strong perturbations: >=25 cells and >=50 AD-BH q<0.05 genes. No energy
  p/q criterion. Features are the union of each strong target's top10 AD
  statistic genes and mean raw UMI >0.25 genes with variance in the dataset's
  top30%. Use mean relative-z response profiles and exact target-feature
  masking to zero where mapping is unambiguous. HDBSCAN correlation distance,
  min_cluster_size4, min_samples1, selection `eom`. Hierarchical leaf ordering
  is for display only, never substituted for the HDBSCAN labels.

## Explicit cross-assay/data limitations

- These are CRISPRko data: the reference's CRISPRi target-RNA repression >=30%
  rule is not an editing-success filter here. DNA indel/frameshift rates remain
  unassessed. Custom ORF targets are not mapped to guessed host genes.
- Reused RNA-QC thresholds are dataset-specific, not the paper's numerical
  mitochondrial/UMI cutoffs. Only available chip/batch identity is used; do
  not invent additional GEM groups or independent biological replicates.
- Singleton NTC guides are not eligible for pairwise core selection because
  a pair of singletons lacks the standardized-AD variance domain. Their
  source rows and exclusion reason remain recorded.
- Positive affine z transformations preserve AD/MW ranks only for a finite
  positive control SD. DE therefore uses UMI-normalized ranks on that defined
  domain; zero-SD genes cannot bypass the reference z representation. The
  first all-feature output remains as a superseded audit, not final DE.
- Agreement with R's asymptotic tail curve establishes numerical equivalence,
  not empirical null calibration for rare/tied genes or highly unbalanced
  samples. Do not promote extreme sparse-gene p values to validated biology.
- The paper's clustering-variance sentence lacks a normalization qualifier;
  raw-UMI variance is the explicit interpretation here, consistent with its
  UMI-mean criterion. This choice is not concealed as verified original code.
- AD calibration lower tail below the grid's near-unity boundary is set to
  p=1 conservatively; extreme machine-underflow tail is bounded at1e-300.
  Off-grid and statistic equivalence tests must precede a full run.
- Local/server binary versions may differ. The execution-host environment is
  recorded and validated; no shared-base upgrade or GPU is authorized.

The 2026 paper is a separate comparison. Its published full methods and
notebook call parameters are not claimed verified from a preprint or function
defaults. Do not quietly mix its pseudobulk/NTC-bootstrap/Leiden paths into
this primary contract. Empty strong-target admission must be reported as
insufficient eligible profiles, not evidence that functional structure is
biologically absent.
