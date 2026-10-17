# Portable NHANES GLP-1 Bridge reconstruction

Version v0.1.0-rc1 is a local offline candidate. It reconstructs the verified descriptive study from eleven manually acquired public NHANES components. No raw XPT, participant-level RDS, historical audit report or precomputed scientific result is supplied as a calculation input.

## Setup

Use the existing software versions in docs/ENVIRONMENT.md. Obtain DEMO_L, HIQ_L, BMX_L, DIQ_L, MCQ_L, KIQ_U_L, BPQ_L, BPXO_L, GHB_L, BIOPRO_L and GLU_L for NHANES August 2021-August 2023 from the official component documentation/download links listed in schemas/public_components.csv. Save each as its documented .xpt filename in a read-only input directory. The candidate verifies locally locked input hashes and labels/units; a changed file fails clearly. It never downloads or writes to inputs.

Copy config/example.json and set PROJECT_ROOT (existing independent workspace), RAW_DATA_DIR, OUTPUT_DIR (explicit new/empty directory outside inputs and code), R_LIBRARY_DIR (optional existing user library), and R_SCRIPT (existing executable). RAW_DATA_DIR/OUTPUT_DIR/library may be absolute or relative to PROJECT_ROOT. Do not leave OUTPUT_DIR empty. Paths with spaces, Unicode or apostrophes are supported by structured argument vectors; use proper shell quoting for the configuration path. Paths involving symlinks/junctions are rejected. Existing nonempty output directories are refused.

Run from a terminal with the configured Python interpreter:

    python -B run.py --config your_config.json

Optionally request display-only figures:

    python -B run.py --config your_config.json --figures

The single entry point verifies candidate hashes and public inputs, then invokes separately auditable modules. No installer, downloader, global library change, unpublished historical RDS or prior manuscript is required. The analytical RDS is generated inside OUTPUT_DIR and must not be redistributed with this code candidate. Inputs are hashed before and after every whole run.

## Execution order

1. Input preflight checks all eleven components, required fields/labels, valid unique identifiers, numeric formats, weights and DEMO alignment.
2. b1/01 builds the analytical object; b1/02,03,04 faithfully calculate archived baseline primary, fasting and conditional results.
3. membership_bridge reuses selected verified baseline assignments to independently calculate membership and pathway counts, replacing private historical membership ledgers.
4. b2/01-07 calculate full-parent survey-domain estimation, D0-D3 routing, exact denominator/common-population tables, age/eGFR, laboratory restrictions, the two paired contrasts, and subgroup uncertainty/variance diagnostics.
5. a3/01-04 verify membership, generate weighted Table1 and raw completeness, and independently check Taylor PSU variance, ratios, interval limits and count identities. Universe denominators come from the reconstructed cohort; expected manuscript counts are not calculation inputs.
6. registry_adapter maps newly calculated CSVs to preserved C02/C03A result identifiers. Metadata mapping contains no estimates. It generates registry/estimands.csv, Table1/completeness and full-precision source tables S1-S7 plus figure source rows.
7. Optional R graphics consume only fresh aggregate panel sources. No new statistical test occurs during display export.
8. tests/concordance.py reads isolated aggregate fixtures only after calculations, compares full precision under frozen tolerances and checks non-estimability/paired identities. A mismatch yields nonzero exit status and is never repaired by changing definitions or widening tolerance.

## Outputs and verification

OUTPUT_DIR contains b1/, b2/, a3/, registry/, source_tables/, presentation/ when requested, logs/, qa/ and provenance/. Logs preserve exact executed commands, module hashes, conditions, R sessions, Python environment, exit codes and input hashes. The default numerical comparisons use exact counts, estimate/SE tolerance1e-8, CI tolerance1e-6, P tolerance1e-12, as documented in tests/fixtures/tolerances.json before test execution. Some unavailable limits are expected NA; they are compared as unavailable rather than replaced by zero or parent-df CIs.

Expected reference fixtures are aggregate-only verification data. They cannot determine a calculated estimate or adapter value and are never imported by R scientific modules or registry_adapter. The original clinical thresholds, missingness, D0/D2 handling, weights and variance policy are retained. Statistical inputs consist solely of public XPTs and freshly generated intermediates.

Read docs/KNOWN_LIMITATIONS.md and docs/LICENSE_DECISION.md. A successful local reproduction does not constitute author approval, public release, actual journal reviewer accessibility, repository archiving or submission readiness.
