#!/usr/bin/env python3
# scCODA on inhibitory-neuron subtype composition (Reviewer 2). Run in sccoda-env.
import os, random, numpy as np, pandas as pd, anndata as ad, tensorflow as tf
from pathlib import Path
from sccoda.util import cell_composition_data as coda_data, comp_ana
SEED=42; os.environ["PYTHONHASHSEED"]=str(SEED); random.seed(SEED); np.random.seed(SEED); tf.random.set_seed(SEED)
BASE=Path("results/abundance"); OUT=BASE/"sccoda"; OUT.mkdir(parents=True,exist_ok=True); SK="orig.ident"
pc=pd.read_csv(BASE/"percell_labels_inh.csv"); cov=pd.read_csv(BASE/"sample_covariates_inh.csv").set_index(SK)
for c in ["genotype","sex","age"]: cov[c]=cov[c].astype(str)
for c in ["med_nCount","med_percent_mt"]: cov[c+"_z"]=(cov[c]-cov[c].mean())/cov[c].std(ddof=0)
F=("C(genotype, Treatment('E33')) + C(sex, Treatment('Male')) + C(age, Treatment('06Mo')) + med_nCount_z + med_percent_mt_z")
obs=pc[[SK,"subtype"]].dropna()
adata=ad.AnnData(X=np.zeros((obs.shape[0],1),dtype="float32"),obs=obs.reset_index(drop=True))
data=coda_data.from_scanpy(adata,cell_type_identifier="subtype",sample_identifier=SK,covariate_df=cov)
m=comp_ana.CompositionalAnalysis(data,formula=F,reference_cell_type="automatic")
r=m.sample_hmc(num_results=8000,num_burnin=2000,num_leapfrog_steps=30,step_size=0.01,verbose=True)
r.effect_df.to_csv(OUT/"sccoda_effects_inh.csv"); r.credible_effects().to_csv(OUT/"sccoda_credible_fdr05_inh.csv")
print("DONE")
