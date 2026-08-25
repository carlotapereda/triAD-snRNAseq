#!/usr/bin/env python3
# Compositional abundance analysis with scCODA (Reviewer 2 compositional concern; Reviewer 3 pt 4:
# model all cell-type counts jointly). Mirrors the lab's sccoda_braak.py pattern.
# Run in sccoda-env:
#   /Users/carlota/miniconda3/envs/sccoda-env/bin/python code/new/abundance_sccoda.py
import os, random
import numpy as np, pandas as pd, anndata as ad, tensorflow as tf
from pathlib import Path
from sccoda.util import cell_composition_data as coda_data
from sccoda.util import comp_ana

SEED = 42
os.environ["PYTHONHASHSEED"] = str(SEED); os.environ["TF_DETERMINISTIC_OPS"] = "1"
random.seed(SEED); np.random.seed(SEED); tf.random.set_seed(SEED)

BASE = Path("results/abundance")
OUT = BASE / "sccoda"; OUT.mkdir(parents=True, exist_ok=True)
SAMPLE_KEY = "orig.ident"

# per-cell labels + per-sample covariates (built by abundance_build_tables.R)
pc  = pd.read_csv(BASE / "percell_labels.csv")
cov = pd.read_csv(BASE / "sample_covariates.csv").set_index(SAMPLE_KEY)
cov["age"] = cov["age"].astype(str)
cov["sex"] = cov["sex"].astype(str)
cov["genotype"] = cov["genotype"].astype(str)
for c in ["med_nCount", "med_percent_mt"]:
    cov[c + "_z"] = (cov[c] - cov[c].mean()) / cov[c].std(ddof=0)

# adjust for genotype, sex, age (treatment-coded) and nuisance covariates (depth, %mt)
FORMULA = ("C(genotype, Treatment('E33')) + C(sex, Treatment('Male')) + "
           "C(age, Treatment('06Mo')) + med_nCount_z + med_percent_mt_z")

def run_level(celltype_col, tag):
    print(f"\n=== scCODA: {tag} ({celltype_col}) ===")
    obs = pc[[SAMPLE_KEY, celltype_col]].dropna().copy()
    obs = obs[~obs[celltype_col].isin(["remove", "negative", "hybrid"])]
    adata = ad.AnnData(X=np.zeros((obs.shape[0], 1), dtype="float32"), obs=obs.reset_index(drop=True))
    data = coda_data.from_scanpy(adata, cell_type_identifier=celltype_col,
                                 sample_identifier=SAMPLE_KEY, covariate_df=cov)
    # 'automatic' reference = the most nearly-constant cell type (avoids arbitrary choice)
    model = comp_ana.CompositionalAnalysis(data, formula=FORMULA, reference_cell_type="automatic")
    res = model.sample_hmc(num_results=8000, num_burnin=2000,
                           num_leapfrog_steps=30, step_size=0.01, verbose=True)
    res.effect_df.to_csv(OUT / f"sccoda_effects_{tag}.csv")
    try:
        res.credible_effects().to_csv(OUT / f"sccoda_credible_fdr05_{tag}.csv")
    except Exception as e:
        print("credible_effects failed:", e)
    print(res.summary())
    print(f"saved sccoda_effects_{tag}.csv")

import sys

LEVELS = {"major": "cell_type_identity", "fine": "hippocampus_cell_type", "cluster": "cluster_id"}

# run every level by default, or only those named on the command line
wanted = sys.argv[1:] or list(LEVELS)
for tag in wanted:
    run_level(LEVELS[tag], tag)
print("\nDONE")
