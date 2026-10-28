# Audit reconstruction only; frozen analyses and assets are never overwritten.
suppressPackageStartupMessages({library(data.table); library(lmerTest)})
root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(root))
out <- file.path(root, "06_reports/2026-10-13_JTM_supplementary_package/provenance")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sp <- fread(file.path(root, "02_intermediate/fig6_visium/spots_galnt_annot.csv.gz"))
sp <- sp[paningrade %in% c("normal", "low_grade", "high_grade") & lib >= 500 & nontissue < 100]
sp[, grade := factor(paningrade, levels=c("normal", "low_grade", "high_grade"))]
sp[, g := as.integer(grade)-1L]
sp[, log10lib := log10(lib)]
arch <- fread(file.path(root, "04_results/2026-10-12_Figure6_GALNT_run01/H1_panin_enzymes_robustness.csv"))
# These two candidates follow the archived spec and the Figure 7 narrative.
forms <- c(acini_collagen="y ~ g + acini + collagen + (1|sample)",
           acini_collagen_depth="y ~ g + acini + collagen + log10lib + (1|sample)")
fits <- list(); checks <- list()
for (i in seq_len(nrow(arch))) {
  gn <- arch$gene[i]; sp[, y := get(gn)]
  for (label in names(forms)) {
    m <- lmer(as.formula(forms[[label]]), data=sp, REML=TRUE)
    cc <- coef(summary(m))["g",]
    a <- arch[i]
    beta_diff <- unname(cc["Estimate"])-a$beta_adj
    p_rel_diff <- abs(unname(cc["Pr(>|t|)"])-a$p_adj)/a$p_adj
    checks[[length(checks)+1L]] <- data.table(gene=gn, model=label, formula=forms[[label]],
      n=nobs(m), n_slides=uniqueN(sp$sample), beta=unname(cc["Estimate"]),
      se=unname(cc["Std. Error"]), df=unname(cc["df"]), p=unname(cc["Pr(>|t|)"]),
      archived_beta=a$beta_adj, archived_p=a$p_adj, beta_abs_difference=abs(beta_diff),
      p_relative_difference=p_rel_diff, matched=abs(beta_diff)<1e-8 && p_rel_diff<1e-5,
      singular=lme4::isSingular(m), REML=TRUE, p_method="Satterthwaite two-sided t test")
    fits[[paste(gn,label,sep="_")]] <- m
  }
}
checks <- rbindlist(checks)
fwrite(checks, file.path(out,"figure6_adjusted_model_comparison.csv"))
saveRDS(fits,file.path(out,"figure6_audit_fits.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"figure6_audit_sessionInfo.txt"))
# Reconstruct all 15 numeric fields without modifying the authoritative CSV.
reconstructed <- rbindlist(lapply(arch$gene, function(gn) {
  sp[, y := get(gn)]
  m <- lmer(y ~ g + (1|sample), data=sp, REML=TRUE)
  ma <- lmer(y ~ g + acini + collagen + log10lib + (1|sample), data=sp, REML=TRUE)
  mc <- lmer(y ~ grade + (1|sample), data=sp, REML=TRUE)
  cc <- coef(summary(m)); ca <- coef(summary(ma)); ct <- coef(summary(mc))
  sm <- sp[,.(v=mean(y)),by=.(sample,grade)]
  sm <- dcast(sm,sample~grade,value.var="v")
  hg <- sm[!is.na(normal) & !is.na(high_grade)]
  lg <- sm[!is.na(normal) & !is.na(low_grade)]
  data.table(gene=gn,beta=cc["g","Estimate"],p=cc["g","Pr(>|t|)"],
    beta_adj=ca["g","Estimate"],p_adj=ca["g","Pr(>|t|)"],
    HG_vs_N=ct["gradehigh_grade","Estimate"],p_HG=ct["gradehigh_grade","Pr(>|t|)"],
    LG_vs_N=ct["gradelow_grade","Estimate"],p_LG=ct["gradelow_grade","Pr(>|t|)"],
    slides_HG_gt_N=sum(hg$high_grade>hg$normal),slides_HG_N=nrow(hg),
    slides_LG_gt_N=sum(lg$low_grade>lg$normal),slides_LG_N=nrow(lg),
    detect_normal=sp[grade=="normal",mean(y>0)],detect_HG=sp[grade=="high_grade",mean(y>0)])
}))
fwrite(reconstructed,file.path(out,"figure6_reconstructed_robustness.csv"))
cell_checks <- rbindlist(lapply(setdiff(names(arch),"gene"),function(v) {
  data.table(gene=arch$gene,field=v,archived=arch[[v]],reconstructed=reconstructed[[v]],
    abs_difference=abs(arch[[v]]-reconstructed[[v]]),
    relative_difference=abs(arch[[v]]-reconstructed[[v]])/pmax(abs(arch[[v]]),1e-300),
    matched=abs(arch[[v]]-reconstructed[[v]]) <= 1e-10*pmax(abs(arch[[v]]),1e-10))
}))
fwrite(cell_checks,file.path(out,"figure6_all_fields_comparison.csv"))
stopifnot(all(cell_checks$matched))
print(checks[,.(gene,model,beta_abs_difference,p_relative_difference,matched,singular)])
cat("Retained spots:",nrow(sp),"slides:",uniqueN(sp$sample),"\n")
cat("Archived numeric fields matched:",nrow(cell_checks),"of",nrow(cell_checks),"\n")
