# Separate guide counting from target-expression annotation

Cell Ranger 8.0.1 accepts the six guide-capture fields (`id`, `name`, `read`, `pattern`, `sequence`, `feature_type`) without the optional `target_gene_id` and `target_gene_name` fields. When those optional fields are supplied, Cell Ranger requires their IDs/names to match the RNA reference. A target absent from that reference must never be renamed `Non-Targeting` to satisfy this check.

`scripts/prepare_guide_counting_reference.py` prepares a new six-column counting reference while preserving the complete original CSV byte-for-byte and a separate authoritative target map. Every guide, its ordering, sequence, capture pattern, read and feature type remain identical. No guide or target is dropped or relabeled. The separately verified downstream control catalogue remains authoritative for NTC/intergenic/positive-control distinctions.

This is suitable for a documented guide-counting/RNA-response workflow that does not request Cell Ranger's optional target-expression annotations. It does not restore missing RNA gene definitions, establish equivalent transcriptome counting between references, or prove guide calls unchanged without execution checks. Downstream preparation must use the preserved original target map, not infer control/target labels from the six-column file. Existing source references and accepted results are never edited.

Run from a clean, pushed revision with absolute source/output paths and the accepted source SHA256:

    python scripts/prepare_guide_counting_reference.py --source /delivery/features.csv --source-sha256 SHA256 --output /private/new-guide-reference --allowed-parent /private --expected-commit COMMIT

The output directory must not exist. `receipt.json` binds the source, unchanged capture fields, preserved target authority, code commit and all payload hashes. Failure does not authorize overwriting a prior attempt or restarting an analysis.
