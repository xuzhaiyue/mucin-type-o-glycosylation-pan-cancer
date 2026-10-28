"""Fig 8 step 2: H8a cluster fraction vs GIScore (spec figure8_spec_v1.md)."""
from __future__ import annotations
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

R = "${PROJECT_ROOT}"
RES = f"{R}/04_results/2026-10-12_Figure8_pathology_run01"


def main() -> None:
    d = pd.read_csv(f"{RES}/H8_slide_table.csv")
    cs = [c for c in d.columns if c.startswith("C") and c[1:].isdigit()]
    d["gi"] = (d.risk - d.risk.mean()) / d.risk.std()
    d["grade"] = d.gradeg.fillna("G1-2")
    rows = []
    for c in cs:
        p = d[c].clip(lower=1e-3)
        d["y"] = np.log(p) - np.log(d[cs].clip(lower=1e-3)).mean(1)  # CLR
        rho, prho = spearmanr(d.risk, d[c])
        m0 = smf.ols("y ~ gi", d).fit(cov_type="HC3")
        m1 = smf.ols("y ~ gi + C(grade) + C(tss_grp)", d).fit(cov_type="HC3")
        m2 = smf.ols("y ~ gi + C(grade)", d).fit(cov_type="HC3")
        hi, lo = d[d.group == "High"][c], d[d.group == "Low"][c]
        rows.append(dict(cluster=c, mean_frac=d[c].mean(), rho=rho, p_rho=prho, beta_unadj=m0.params["gi"], p_unadj=m0.pvalues["gi"],
                         beta_grade=m2.params["gi"], p_grade=m2.pvalues["gi"],
                         beta_adj=m1.params["gi"], se_adj=m1.bse["gi"], p_adj=m1.pvalues["gi"], lo=m1.conf_int().loc["gi", 0], hi=m1.conf_int().loc["gi", 1],
                         p_MWU_high_low=mannwhitneyu(hi, lo).pvalue, mean_high=hi.mean(), mean_low=lo.mean()))
    t = pd.DataFrame(rows)
    for k in ("p_rho", "p_unadj", "p_grade", "p_adj"):
        t["bh_" + k[2:]] = multipletests(t[k], method="fdr_bh")[1]
    t.to_csv(f"{RES}/H8a_cluster_vs_GIScore.csv", index=False)
    pd.set_option("display.width", 250); print(t.round(4).to_string())
    # TSS composition of clusters (confound view)
    ct = d.groupby("tss_grp")[cs].mean(); ct["n"] = d.groupby("tss_grp").size(); ct.to_csv(f"{RES}/H8c_cluster_by_TSS.csv")
    print(ct.round(3))
    # how much variance is TSS vs GIScore (R2 of TSS only on each cluster CLR)
    out = []
    for c in cs:
        d["y"] = np.log(d[c].clip(lower=1e-3)) - np.log(d[cs].clip(lower=1e-3)).mean(1)
        out.append(dict(cluster=c, r2_tss=smf.ols("y ~ C(tss_grp)", d).fit().rsquared, r2_gi=smf.ols("y ~ gi", d).fit().rsquared))
    pd.DataFrame(out).to_csv(f"{RES}/H8c_r2_tss_vs_giscore.csv", index=False); print(pd.DataFrame(out).round(3))
    print("GIScore by TSS ANOVA R2:", smf.ols("risk ~ C(tss_grp)", d).fit().rsquared)


if __name__ == "__main__":
    main()
