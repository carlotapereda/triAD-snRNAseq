#!/usr/bin/env python3
# Beta-GLM on inhibitory-subtype proportions with factorial + nuisance. Run in pymc-env.
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pymc as pm
from pathlib import Path
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
BASE=Path("results/abundance"); OUT=BASE/"glm"; OUT.mkdir(parents=True,exist_ok=True); SK="orig.ident"
cov=pd.read_csv(BASE/"sample_covariates_inh.csv").set_index(SK)
cov["nCount_z"]=zscore(cov["med_nCount"]); cov["mt_z"]=zscore(cov["med_percent_mt"])
cov["genoE44"]=(cov["genotype"]=="E44").astype(float); cov["sexF"]=(cov["sex"]=="Female").astype(float)
cov["age12"]=(cov["age"]=="12Mo").astype(float); cov["age18"]=(cov["age"]=="18Mo").astype(float)
cov["geno_sex"]=cov["genoE44"]*cov["sexF"]; cov["geno_age12"]=cov["genoE44"]*cov["age12"]; cov["geno_age18"]=cov["genoE44"]*cov["age18"]
TERMS=["Intercept","genoE44","sexF","age12","age18","geno_sex","geno_age12","geno_age18","nCount_z","mt_z"]
REPORT=["genoE44","sexF","age12","age18","geno_sex","geno_age12","geno_age18"]
X=np.column_stack([np.ones(len(cov))]+[cov[t].values for t in TERMS[1:]]).astype(float)
cnt=pd.read_csv(BASE/"counts_inh.csv").set_index(SK).loc[cov.index]; props=cnt.div(cnt.sum(axis=1),axis=0)
def fit(y):
    with pm.Model():
        b=pm.Normal("beta",0,1.5,shape=X.shape[1]); phi=pm.HalfNormal("phi",25.); mu=pm.math.sigmoid(pm.math.dot(X,b))
        pm.Beta("obs",alpha=mu*phi,beta=(1-mu)*phi,observed=np.clip(y,1e-5,1-1e-5))
        tr=pm.sample(draws=2000,tune=2000,chains=4,target_accept=0.95,random_seed=42,progressbar=False,cores=2)
    post=tr.posterior["beta"].stack(s=("chain","draw")).values; out=[]
    for n in REPORT:
        s=post[TERMS.index(n),:]; p=2*min((s>0).mean(),(s<0).mean())
        out.append(dict(term=n,coef=float(s.mean()),p_val=float(p)))
    return out
rows=[]
for ct in props.columns:
    print("->",ct)
    for r in fit(props[ct].values): rows.append(dict(subtype=ct,**r))
df=pd.DataFrame(rows); df["p_adj"]=multipletests(df["p_val"],method="fdr_bh")[1]; df["sig"]=df["p_adj"]<0.05
df.to_csv(OUT/"glm_inh.csv",index=False)
print(df[df["sig"]].sort_values("p_adj").to_string(index=False) or "no sig"); print("DONE")
