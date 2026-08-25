# Tri-AD: snRNA-seq of the aging APOE knock-in mouse hippocampus

Analysis code for Li et al., "Tri-AD" (Nature Communications, in revision). Single-nucleus
RNA-seq of hippocampus from humanized APOE3/3 and APOE4/4 knock-in mice, male and female, at
6, 12, and 18 months (12 conditions, 47 libraries, 555,008 nuclei x 27,153 genes).

Processed data can be explored and downloaded at https://ucsftriad.org.
Archived release: [ZENODO DOI].

## Software versions

### Original pipeline (scripts 1 to 6, run January 2023)

- Cell Ranger 4.0.0 (`cellranger count`, default parameters) against a custom mm10-1.2.0
  reference including introns. Alignment is not repeated by these scripts; they start from the
  per-library count matrices.
- Seurat v4.0.2 (version recorded inside the deposited Seurat object), sctransform v2
  (`SCTransform(vst.flavor = "v2")`).
- The full command log recorded in the deposited object (function calls, parameters, seeds,
  timestamps) is provided in `object_provenance.txt`. Key parameters: PCA computed with the
  default `npcs = 50`; UMAP, neighbor graph, and clustering used `dims = 1:15`; Louvain
  clustering at `resolution = 0.5`, giving 36 clusters.

### Revision analyses (added during peer review)

- R 4.5.2 with Seurat 5.5.0. Exact versions of every R package are pinned in `renv.lock`;
  recreate the environment with `renv::restore()` (install renv first). The full
  `sessionInfo()` output is in `sessionInfo.txt`.
- MAST 1.36.0 (Bioconductor, R 4.5.3) in its own conda environment: `envs/mast-env.yml`.
- Python abundance models in two conda environments: `envs/sccoda-env.yml` (scCODA 0.1.9) and
  `envs/pymc-env.yml` (PyMC 5.26).

```bash
conda env create -f envs/sccoda-env.yml
conda env create -f envs/pymc-env.yml
conda env create -f envs/mast-env.yml
```

Note: the merged object is large. Loading it requires roughly 20 GB of RAM; on machines with
24 GB, raise the R vector memory limit at runtime with `mem.maxVSize(vsize = 60000)`.

## Workflow and scripts

Run the numbered scripts in order. Each script carries a header comment stating its purpose,
inputs, and outputs.

| # | Script | Purpose | Inputs | Outputs |
|---|--------|---------|--------|---------|
| 1 | `01_create_obj_and_QC.R` | Read per-library counts, merge, QC-filter, SCTransform (v2), PCA, UMAP, Louvain clustering (res 0.5, 36 clusters) | 47 Cell Ranger count matrix directories (`Read10X`) | Merged, filtered, clustered Seurat object (`Emouse.rds`); QC violin plots |
| 2 | `02_cell_type_annotation.R` | Score marker modules (`AddModuleScore`) and assign per-nucleus cell type labels; flag hybrid/negative nuclei; assign cluster identities from majority label | `Emouse.rds` | Metadata columns `cell_type_gen`, `Seurat_Clusters` |
| 3 | `03_DEG_analysis.R` | APOE4/4 vs APOE3/3 differential expression per cell type, stratified by sex and age (`FindMarkers`, MAST) | Annotated object | All-gene DEG table (`DEG.list.all`) |
| 4 | `04_cell_type_abundance.R` | Per-library cell type proportions (hybrid/negative removed); three-way Type III ANOVA (genotype x sex x age) | Annotated object | Proportion tables; `2025oct_celltype_3anova_results.csv` |
| 5 | `05_cluster_sample_abundance.R` | Same design at the level of the 36 clusters | Annotated object | Proportion tables; `cluster_3anova_results.csv` |
| 6 | `06_APOE_expression_level_comparison.R` | hApoE transgene expression summaries (SCT residuals) per genotype-sex-age group; pairwise tests | Annotated object | `APOEexpr_pairwise_pval.csv` |
| — | `requested_code.R` | Regenerate DEG tables (`DEG.csv`, `sDEG.csv`, `Inh_all_gene_SCT_DEGs.csv`), UpSet plots, and pathway heatmaps | Annotated object; pathway CSVs | DEG tables and figure panels |

### QC filters (script 1)

Nuclei were retained with more than 250 detected genes, at least 500 counts, and
log10GenesPerUMI above 0.85; genes expressed in fewer than 10 nuclei were removed. These
thresholds were verified against the QC metric minima of the deposited object (555,008 nuclei).

### Revision scripts

Analyses added during peer review, in `revision/`:

| Script | Purpose | Inputs | Outputs |
|--------|---------|--------|---------|
| `QC_validation.Rmd` | Sample mixing (iLISI/cLISI), QC distributions, annotation marker panels | `Emouse.rds` | QC and mixing figures |
| `resolution_sweep.R`, `stability_fig_from_labels.R` | Reclustering at resolutions 0.3 to 1.0; ARI and nestedness | `Emouse.rds` | Stability figures |
| `umap_and_elbow.R` | Labeled UMAP, vascular highlight, PCA elbow | `Emouse.rds` | UMAP and elbow figures |
| `subtype_and_oligo_validation.R` | Oligodendrocyte identity, DAM/DAA/DAO scores, microglia and vascular subsets | `Emouse.rds` | Validation figures and tables |
| `export_metadata.R` | Export per-nucleus metadata so downstream steps avoid reloading the 19 GB object | `Emouse.rds` | `cell_metadata.rds` |
| `abundance_build_tables.R`, `inh_build_tables.R` | Per-library count tables for abundance models | `cell_metadata.rds` | Count tables (CSV) |
| `abundance_sccoda.py`, `inh_sccoda.py` | scCODA compositional abundance models (nuisance-adjusted) | Count tables | Credible-effect tables, forest plot inputs |
| `abundance_glm.py`, `inh_glm.py` | Beta-binomial regression with full factorial design | Count tables | Model summaries |
| `object_provenance.R` | Extract the command log and versions recorded in the object | `Emouse.rds` | `object_provenance.txt` |

## Notes on reproducibility

- Random seeds are the Seurat defaults and are recorded, with all other parameters, in
  `object_provenance.txt`. UMAP coordinates may differ cosmetically across platforms; cluster
  memberships do not.
- The marker module named "Endothelial" in script 2 scores Csf1r, a myeloid marker; this module
  therefore captures myeloid rather than endothelial identity. The code is deposited as run;
  see the response to reviewers for the validation analyses that address this.
- `requested_code.R` contains four corrections relative to the version used originally
  (documented in its header): the sDEG filtering step, a column name, a cell type label, and a
  heatmap NA filter. The corrected NA filter reproduces the published Supp Fig 4 term set
  exactly.

## Citation

[Full citation once available.]
