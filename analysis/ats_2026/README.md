# ATS 2026 central analyses and figures

This directory preserves the pooled ASE Results paragraph and final figure scripts from October 22, 2026. These are central aggregate analyses; sites do not run them as part of the site pipeline.

The ASE analysis replaces Rush's earlier contribution with completed run `20261022_204133_4103978` and retains Hopkins and the other systems once. The cohort includes 701,028 hospitalizations, 132,111 ASE and 56,215 presumed-infection-without-ASE hospitalizations. The original 13-organism and 77-pair comparison families are retained. The clinical figure displays 29 selected organism–antimicrobial pairs; that selected figure family is separate from the 77-pair inference family. Organism models exclude Emory because of unresolved result classification concerns.

## Reproduce with the prepared aggregate workspace

Set `CLIF_ANALYSIS_ROOT` to the local aggregate analysis workspace. It must contain the prepared directories `ASE_RUSH_updated_20261022/ase/inputs`, `ASE_RUSH_updated_20261022/clinical_figure`, `ASE_RUSH_updated_20261022/monthly_main_analysis_input.csv`, and `CLIF_refresh_20261022/atlas`. The ASE directory also needs `original_77_pair_screen.csv` and `susceptibility_screen_functions.R`. The workspace source registries and validation records identify the exact input snapshots. Aggregate input exports are deliberately not copied into Git; these scripts require that prepared workspace and do not fetch Box data.

Run from the repository root, in order:

```sh
Rscript --vanilla analysis/ats_2026/ase_analysis.R
Rscript --vanilla analysis/ats_2026/ase_standardized_profiles.R
Rscript --vanilla analysis/ats_2026/ase_clinical_pairs.R
Rscript --vanilla analysis/ats_2026/ase_figure.R
Rscript --vanilla analysis/ats_2026/atlas_figure.R
```

Outputs are written to the prepared aggregate workspace. The checked-in ASE PNG/PDF and atlas vector PDF are the reviewed submission figures. The atlas script also exports the 24,000-pixel-wide raster in that workspace. Figure artifacts in this directory are snapshots; copy regenerated figures here explicitly after review.

## Interpretation and outstanding QC

The display labels Non-infected and Infected mean not meeting presumed-infection criteria and presumed infection without ASE, respectively. These are retrospective hospitalization classifications, not independent clinical infection diagnoses. Resistance is the requested display terminology; the measured outcome remains CLIF non-susceptibility among interpretable results.

Available uncertainty bars are site-clustered 95% confidence intervals; daggered clinical pairs are descriptive pooled fractions without model-derived intervals. Models standardize across shared observed site-month support, with equal site weights. Missing observations are not assigned zero.

Rush's replacement upload resolves the earlier extremely low classification counts and passes cohort reconciliation. Previously identified fungal–antibacterial and bacterial–antifungal pairings in source susceptibility exports remain unresolved; susceptibility findings remain provisional pending source linkage review. Removing fungi from the atlas resistance wheel is a display decision and does not validate remaining pairs. Fungi remain on the detection wheel. The 47-hospital scope in the Results is based on the supplied consortium inventory.
