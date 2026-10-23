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
  It never grants permissions or accepts a requested mode change. `LD_PRELOAD`
  must be set for the owned run only, not in a shell profile/system setting.
- `stage_view.py` creates an unused native software-only view. It hash-verifies
  independent copies of the original launcher and candidate runtime and
  references unchanged installed components. This is not a large-data workspace.

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

Never enable automatic recovery of a real failed analysis as part of this test.
Do not modify shared mounts, permissions, original FASTQ/BAM, previous outputs,
or the installed vendor package. Record runtime, shim, inputs and outputs with
hashes. Full production recovery needs a separate progress/cost assessment.
