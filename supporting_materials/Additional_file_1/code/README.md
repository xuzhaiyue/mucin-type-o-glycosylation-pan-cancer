# Analysis code accompanying the article

Archived primary scripts are in `archived/`; their original bytes and hashes are preserved. `PIPELINE.md` records their ordering. This is a scientific code snapshot, not a newly retrained model or a fully portable one-command release.

The original scripts include historical absolute paths and write to their original analysis folders. Configure a separate writable project workspace before execution. Review `archived_code_path_inventory.csv`, supply public data and the source annotation objects, and redirect output paths. Do not execute training or plotting scripts in the frozen project.

`../figure6_audit/verify_figure6.R` is a new, explicitly identified audit reconstruction of the Figure 6E robustness table. It consumes the archived per-spot GALNT table, keeps the original 608-spot gate and reproduces all 84 numerical entries. Original scripts 31/32/33 remain available; only the archived robustness-table writer is missing.

The original whole-study environment lock was not recovered. The included sessionInfo is for this audit only. GIScore coefficients/model hashes and all scientific results remain frozen. External ABMIL seed timing and earlier single-cell annotation dependencies limit exact full-pipeline reproducibility.

Publication numbering: historical S4 → S1; S5 → S2; S6 → S3; S8 → S4. Historical source-table names and scripts are retained as provenance.
