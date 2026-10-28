# Public repository release notes

This repository stages source code, derived numerical tables, specifications, model-lock records, figures, and provenance files associated with the submitted manuscript. The original local project remains unchanged.

Machine-specific absolute paths in text files were replaced with environment-style placeholders (for example, `${PROJECT_ROOT}` and `${PANIN_ROOT}`) to avoid exposing workstation paths. Numeric source-table values and figure content were not altered. Package-level `SHA256SUMS.txt` manifests were recomputed after this path-only cleanup. The local preparation report retains the checksums of the exact source-package members before cleanup.

Several archived scripts refer to raw inputs that are not redistributed here and are not a single turnkey workflow. Configure the path placeholders and obtain upstream inputs under their applicable access terms before attempting a rerun. The manuscript reports archived outputs and documents the limits to exact retraining.

The journal workbook contains 55 selected figure-linked tables (5,931 rows) for readability. The complete 123 CSV/CSV.GZ source tables (620,339 rows) remain in the repository. No underlying numerical records were removed from the public release.

No software or data license was supplied for this release; upstream dataset and provider terms continue to apply.
