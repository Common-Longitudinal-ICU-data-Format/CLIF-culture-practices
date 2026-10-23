# Buddy testing before multisite use

## Run and execution checks

Use a site with both `microbiology_culture` and `microbiology_susceptibility`. Copy `config/config_template.json` to the ignored `config/config.json`, and set the site label, local repository/table paths, and file type. The study team maintains the shared dates in `utils/study_settings.R` (currently January 1, 2018–December 31, 2024); culture coverage is assumed complete within that window.

Before the first test or pipeline run, install the pinned packages from the repository root:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
```

This is a separate setup step. The pipeline does not install missing packages. The project environment activates automatically thereafter; restore again when `renv.lock` changes or preflight reports a dependency mismatch.

Then run the checks and pipeline:

```sh
Rscript tests/test_preflight.R
Rscript tests/test_config.R
Rscript tests/test_core.R
Rscript tests/test_pooling.R
Rscript tests/run_integration.R
Rscript code/00_run_pipeline.R
```

Check the completed run manifest and empty privacy audit. Keep private intermediates local. Review `susceptibility_availability`, `susceptibility_qc`, `monthly_susceptibility_linkage_qc`, monthly organism–drug aggregates, temporal models, fitted curves, and coverage figures from that same run. Synthetic fixtures verify execution and counting rules; clinical source reconciliation is still required. UCMC has no AST table and cannot perform that reconciliation.

## Reconcile source laboratory records

1. Compare culture events and positive isolates by month and specimen against the source laboratory system. Review sudden disappearances and category changes; an organism absence is interpretable only with complete culture capture and stable identification/mappings. Verify that simultaneous specimens are not collapsed incorrectly by the event key.
2. Trace a sample of `organism_id` links in private intermediates. Include polymicrobial specimens, repeated records, missing IDs, duplicate IDs across events, and unmatched AST records. Missing-ID isolates must remain in coverage denominators. A two-linked/eight-unlinked sample should show 20% linkage and 20% testing when both linked isolates are tested, even though linkable-only testing is 100%.
3. Reconcile canonical organism, antimicrobial and S/NS categories with the pinned mCIDE. Verify susceptible and non-susceptible reports separately; raw MIC/text never supplies a missing standardized interpretation. Confirm how the local ETL maps intermediate/SDD and breakpoint revisions. Conflicting S and NS records for one isolate/drug become indeterminate; the schema has no timestamp to adjudicate them.
4. Compare S, NS, indeterminate, unavailable and unreported-test counts for selected organism–drug–specimen months with source reports. The tested-fraction denominator is S + NS; the testing-coverage denominator is all observed positive isolates of that organism. An organism may be S to one drug and NS to another; no universal resistance label or inferred MRSA/VRE/ESBL/CRE phenotype is assigned by this analysis.
5. Test months with organisms but no interpretable AST, and months with zero organisms while cultures continue. The former must stay unknown. The latter enter detection-rate models under the shared complete-coverage assumption; their tested fraction stays missing. A month without specimen-source activity is excluded.

Report source or mapping gaps to the study team so analysis restrictions can be coordinated. Defaults are 90% linkage completeness and 50% testing coverage for rate models; these are QC screens. Inspect losses at these thresholds and compare plausible stricter thresholds before drawing conclusions.

## Scientific sensitivity checks

These checks require site laboratory context and additional local reruns; they are not yet automated analysis branches.

- **First isolate:** compare primary repeated-isolate detection trends against the first isolate of each organism per patient/hospitalization, using a prespecified specimen policy. This assesses whether frequent reculturing drives a trend. Do not deduplicate independently by drug in a way that selects results.
- **Stable panel:** restrict to organism–drug combinations and intervals with stable testing policies, breakpoints and mappings. Compare coverage and fitted effects with primary outputs. Selective/cascade reporting can change the tested fraction even when overall coverage appears stable.
- **Case/specimen mix:** compare overall and specimen-stratified results and review laboratory/platform or service changes. Neither ICU-day normalization nor seasonal adjustment controls clinical case mix.

Before pooling, register only completed runs with matching code/mCIDE hashes and source-validated dates. Set `culture_qc_pass` and `ast_qc_pass` after review. Preserve the review record locally. Report sparse-model exclusions, coverage changes, FDR and residual diagnostics alongside any increase/decrease label.

Dependency setup uses the committed `.Rprofile` and `renv/activate.R` bootstrap (renv 1.1.5). Restore with `Rscript -e 'renv::restore(prompt = FALSE)'` from the repository root. The site runner checks dependencies and table schemas before analysis, and writes aggregate date availability to `provenance/preflight_source_availability.csv` in the completed run.

## Infection and ASE hospitalization subgroups

The pipeline also runs `code/11_ase_stratified_analysis.R`, using full hospitalization clinical data to classify adult hospitalizations into three mutually exclusive groups: No presumed infection, Presumed infection without ASE, and ASE. Each group has its own ICU admission/day denominators, including uncultured stays, with additional culture, organism, timing, and optional susceptibility outputs. Missing required ASE inputs skip this additional analysis explicitly. Restore the updated `renv.lock` before running; DuckDB/DBI are included and no Python setup or additional site configuration is needed. See [ASE analysis definitions](ase_analysis.md).


## Troubleshooting susceptibility linkage

If script 9 stops with `organism_id maps to multiple ICU isolate events or categories`, the same ID has more than one retained culture-event, patient, hospitalization, organism category, specimen category, or collection timestamp. Matching raw-label variants of one identified isolate are collapsed without multiplying susceptibility counts. Genuine conflicts still stop the analysis; do not resolve them by keeping an arbitrary first row or by joining every conflicting culture to the same susceptibility record.

The failed run writes an aggregate `susceptibility_linkage_conflict_summary` under its `susceptibility/` output folder. That summary reports counts of conflicting IDs by dimension and can be shared with the study team. The detailed `organism_id_conflicts` and `organism_id_conflict_culture_rows` files are saved under `data/intermediate/runs/<run_id>/susceptibility/`; these contain identifiers and exact timestamps and must stay local. Review those records against the source specimen/isolate accession and the ETL key shared by both microbiology tables. Repeated collections of the same species need distinct isolate identifiers; susceptibility rows for different drugs tested on one isolate should share its identifier. Distinct order/event records for one underlying specimen also require source reconciliation before assigning results.
