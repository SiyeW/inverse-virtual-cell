# Experimental Cell Ranger 8.0.1 CIFS compatibility

Not a production-ready analysis release. This adapter is opt-in and does not
change aligners, cell calling, UMI counting, or guide assignment. Existing
installations and scientific attempts must remain unchanged.

## Components

- Upstream Martian v4.0.12, commit `2826c42e541ad2929cbb9ef73fb84c9fc341943f`.
  Its built-in `MRO_UNIQUIFIED_DIRECTORIES=disable` uses the existing legacy
  directory path. It is not a new implementation of attempt/resume handling.
- `martian-v4.0.12-cifs.patch` adds **only** opt-in final publication using
  regular directories and same-filesystem hard links. Original files remain
  in place; collisions require matching content; unsupported/special/external
  output paths fail. Enable with `MRO_CIFS_COMPAT=1`.
- `same_mode_fchmod.c` handles an observed Rust `fs::copy` failure on CIFS.
  It first calls the real operation. Only EPERM/EACCES, an exact match between
  current and requested permission bits, a CIFS descriptor, and a path beneath
  `MRO_CIFS_SAME_MODE_ROOT` permit the already-satisfied operation to succeed.
  It never grants permissions or accepts a requested mode change. The original
  launcher clears inherited `LD_PRELOAD`. The isolated view sets its private
  `.env.json` to the hash-verified shim inside that view, not to a caller-selected
  library. No shell profile or system loader configuration is changed.
- `stage_view.py` creates an unused native software-only view. It hash-verifies
  independent copies of the original launcher and candidate runtime and
  references unchanged installed components. This is not a large-data workspace.
- `copy_fallback.py` restricts only the owned process tree: the Linux x86_64
  `copy_file_range` syscall returns ENOSYS, selecting the existing ordinary
  read/write fallback. This addresses an observed CIFS EIO from that optional
  optimization. All other syscalls are unchanged; no system policy or mount is
  modified, and privileges cannot be gained in the child tree.

## Log identity recovery

Upstream `VerifyLogFile` compares the open descriptor's file identity with the
current `_log` path. CIFS can invalidate that comparison without deleting the
log. `patch_log_identity.py` is a SHA-gated adapter for the pinned logging source;
apply it only to a fresh software-source copy, then build a separate runtime.

With `MRO_CIFS_COMPAT=1`, an explicitly bounded `MRO_CIFS_SAME_MODE_ROOT`, and a
kernel-confirmed CIFS descriptor, an identity mismatch can use a live-content
proof. Under the logger mutex, it appends an unpredictable marker through the
original descriptor, syncs, and reads that marker at its actual append offset
through the expected path. Missing files, symlinks, scope escapes, filesystem
changes, and actual replacement files fail closed, including a replacement with
an identical old prefix. Native filesystems and ordinary runs retain the strict
upstream check. This changes log checking, not scientific stages or their inputs.

Copy `cifs_log_identity*.go` into the candidate's `martian/util` through the
adapter. Run util/core/mrp tests and the explicit actual-storage case via
`MRO_CIFS_LOG_TEST_ROOT`. The remapped-inode case simulates an unequal ID; it is
not a claim to reproduce the original multi-day kernel event. The actual CIFS
case must also reject a real copied replacement. Re-run the small interrupted
pipeline fixture with the new runtime before a separately authorized recovery.

`resume_preserved.py` is a one-attempt production recovery entrypoint. It requires
an immutable stopped-run inventory, unchanged invocation/UUID, exact scientific
CSV and input identities, a tested private view, and explicit authorization to
preserve completed chip02 work. Completed metadata are hashed; existing outputs
and required remaining-stage inputs are checked by their recorded metadata,
without pretending to freshly hash hundreds of GB. Completed-stage temporary
outputs that have already been consumed are recorded separately and cannot
replace validation of the actual remaining-stage dependencies. Original control
logs remain in place; the resumed owner uses a fresh control directory. It
refuses completed-node execution and upstream shard/barcode reexecution, keeps
the existing CPU/RAM/disk protections, and never enables full-stage reset.

## Build and tests

Use an independent checkout of the fixed upstream commit. Export with LF line
endings (`git -c core.autocrlf=false archive` on Windows), apply the patch, and
build with a recorded Go toolchain. Keep upstream license `MARTIAN_LICENSE`.
Dependencies must match upstream `go.sum`; offline mirrors are transport only,
not an authority replacing the pinned checksums.

Run upstream core tests, `TestCifs*` on both native and actual CIFS storage,
`test_stage_view.py`, and the explicit `test_same_mode.py` on a fresh disposable
CIFS root. The symlink-rejection test requires native symlink support; it may
skip on CIFS but must pass on native storage. CIFS can report different inode
numbers for aliases, so hard-link tests verify shared synthetic contents and
survival after one name is removed, not inode equality alone.

Then run the complete official tiny fixture with no BAM, compare exact count
arrays/barcodes/features with an accepted unmodified native run, and test
interrupt/resume. A version/help check, preflight pass, or unit-test result does
not establish end-to-end or production readiness. For CRISPR use, an additional
bounded RNA+guide fixture must also preserve original guide identities.

`run_resume_fixture.py` is a synthetic interruption/resume check. Its first
stage can execute only once; the second stage is deliberately interrupted.
Successful recovery must reuse the completed step and publish the value 42.
It accepts only a new provenance test directory, never an existing analysis.
`verify_fixture.py` also rejects publication errors and malformed final metadata
even when the vendor program returns exit code zero.

`replay_prefix_fixture.py` admits only the previously accepted 100,000-pair
RNA/CRISPR prefixes and engineering-only forced 10-cell recipe. It verifies
those small derived FASTQ hashes, the reference's frozen metadata, baseline
output hashes, and the independent software view before execution. It never
regenerates prefixes or reads full original FASTQs. Success additionally requires
exact matrix arrays, guide calls, feature definitions, QC and output coverage.

The optional report overlay changes only the attachment copy policy: verified
contents inherit destination-filesystem permissions rather than copying a local
mode onto CIFS. `adapt_report.py` is a narrow, source-hash-gated transformer.
Vendor source and its generated derivative stay in the internal installation;
they are not redistributed here. Cell Ranger's own internal-use license is
distinct from Martian's MIT license. The original vendor files are preserved.

The explicit `--resume-report-only` option is limited to the diagnosed failed
100k-prefix reporting step. It archives prior evidence with hashes, requires the
same invocation and stops if counting/guide assignment unexpectedly runs again.
It is not a general recovery entrypoint for production runs.

Never enable automatic recovery of a real failed analysis as part of this test.
Do not modify shared mounts, permissions, original FASTQ/BAM, previous outputs,
or the installed vendor package. Record runtime, shim, inputs and outputs with
hashes. Full production recovery needs a separate progress/cost assessment.
