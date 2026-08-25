#!/usr/bin/env python3
# Beta-GLM for cell-type proportions with the full age x sex x genotype factorial (main effects AND
# interaction terms) plus nuisance covariates (sequencing depth, %mito). Addresses Reviewer 3 pt 3
# (joint GLM of primary + nuisance covariates with interaction terms, at the library level so there is
# no cell-level pseudo-replication). Mirrors the lab's PyMC beta-model pattern.
# Run in pymc-env:
#   /Users/carlota/miniconda3/envs/pymc-env/bin/python code/new/abundance_glm.py
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pymc as pm
from pathlib import Path
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests

BASE = Path("results/abundance"); OUT = BASE / "glm"; OUT.mkdir(parents=True, exist_ok=True)
SAMPLE_KEY = "orig.ident"
N_DRAWS, N_TUNE, TARGET_ACCEPT, CHAINS, SEED = 2000, 2000, 0.95, 4, 42

cov = pd.read_csv(BASE / "sample_covariates.csv").set_index(SAMPLE_KEY)
cov["nCount_z"] = zscore(cov["med_nCount"]); cov["mt_z"] = zscore(cov["med_percent_mt"])
# design: genotype (E33 ref), sex (Male ref), age (06Mo ref), + interactions + nuisance
cov["genoE44"] = (cov["genotype"] == "E44").astype(float)
cov["sexF"]    = (cov["sex"] == "Female").astype(float)
cov["age12"]   = (cov["age"] == "12Mo").astype(float)
cov["age18"]   = (cov["age"] == "18Mo").astype(float)
cov["geno_sex"]   = cov["genoE44"] * cov["sexF"]
cov["geno_age12"] = cov["genoE44"] * cov["age12"]
cov["geno_age18"] = cov["genoE44"] * cov["age18"]

TERMS = ["Intercept","genoE44","sexF","age12","age18","geno_sex","geno_age12","geno_age18","nCount_z","mt_z"]
REPORT = ["genoE44","sexF","age12","age18","geno_sex","geno_age12","geno_age18"]  # skip nuisance in FDR table

def design(df):
    X = np.column_stack([np.ones(len(df))] + [df[t].values for t in TERMS[1:]]).astype(float)
    return X

def clamp01(y, eps=1e-5): return np.clip(y, eps, 1 - eps)

def fit(y, X):
    with pm.Model():
        beta = pm.Normal("beta", 0, 1.5, shape=X.shape[1])
        phi  = pm.HalfNormal("phi", 25.0)
        mu   = pm.math.sigmoid(pm.math.dot(X, beta))
        pm.Beta("obs", alpha=mu*phi, beta=(1-mu)*phi, observed=clamp01(y))
        tr = pm.sample(draws=N_DRAWS, tune=N_TUNE, chains=CHAINS, target_accept=TARGET_ACCEPT,
                       random_seed=SEED, progressbar=False, cores=2)
    post = tr.posterior["beta"].stack(s=("chain","draw")).values  # (nterms, ndraw)
    out = []
    for name in REPORT:
        s = post[TERMS.index(name), :]
        p = 2 * min((s > 0).mean(), (s < 0).mean())
        out.append(dict(term=name, coef=float(s.mean()),
                        hdi2_5=float(np.quantile(s, .025)), hdi97_5=float(np.quantile(s, .975)),
                        p_val=float(p)))
    return out

def run_level(counts_file, tag):
    print(f"\n=== beta-GLM: {tag} ===")
    cnt = pd.read_csv(BASE / counts_file).set_index(SAMPLE_KEY)
    cnt = cnt.loc[cov.index]
    props = cnt.div(cnt.sum(axis=1), axis=0)
    X = design(cov)
    rows = []
    for ct in props.columns:
        print("  ->", ct)
        for r in fit(props[ct].values, X):
            rows.append(dict(level=tag, celltype=ct, **r))
    df = pd.DataFrame(rows)
    df["p_adj"] = multipletests(df["p_val"], method="fdr_bh")[1]
    df["sig"] = df["p_adj"] < 0.05
    df.to_csv(OUT / f"glm_{tag}.csv", index=False)
    print(df[df["sig"]].sort_values("p_adj").to_string() or "  (no significant terms)")
    return df

a = run_level("counts_major.csv", "major")
b = run_level("counts_fine.csv", "fine")
print("\nDONE")
