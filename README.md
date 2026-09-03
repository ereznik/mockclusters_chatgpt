# Mixed-type patient clustering benchmark

This project accepts any patient-by-feature CSV, executes eight mixed-data
clustering strategies, and creates a heatmap-heavy PDF report. It also includes a
400-patient, 100-feature mock cohort with four known clusters, correlated noise,
outliers, and structured missingness.

## Quick start

Start with the real-data-shaped example file:

```text
data/example_patients_input.csv
```

It is a runnable, entirely synthetic example of the exact layout expected for
your own data. It has one header row, one patient per subsequent row, one unique
ID column, and only candidate clustering features after the ID. It deliberately
does **not** contain `true_cluster` or any other answer column.

Run all methods and create the PDF with:

```bash
Rscript analyze_csv.R \
  --input data/example_patients_input.csv \
  --output-dir output/example_patients \
  --id-column patient_id \
  --k auto
```

The PDF is always written to `output/pdf/`. For this file, the default report is
`output/pdf/example_patients_clustering_report.pdf`; a prebuilt copy is committed
there.

To use your own data, save it anywhere as a CSV and change only the paths:

```bash
Rscript analyze_csv.R \
  --input path/to/your_patients.csv \
  --output-dir output/my_patients \
  --id-column patient_id \
  --k auto
```

The command infers a feature dictionary if none is supplied and writes it to the
run directory for review. Its PDF is saved as
`output/pdf/your_patients_clustering_report.pdf`. Use `--pdf a_name.pdf` to choose
a different filename; directory components are ignored so the file still goes
to `output/pdf/`. Use `Rscript analyze_csv.R --help` for every option.

Reproduce the mock-data analysis and final PDF with:

```bash
Rscript run_all.R
```

The only statistical package required by the analysis is `cluster`, which is a
recommended R package. The PDF report additionally uses `ggplot2` and its normal
dependencies.

The final report is written to:

```text
output/pdf/mixed_clustering_benchmark.pdf
```

## Exact CSV structure for your data

Use `data/example_patients_input.csv` as the model. Its first columns and rows
look like this (the file contains 60 rows so it can be run immediately):

```csv
"patient_id","age_years","bmi_kg_m2","systolic_bp_mmhg","crp_mg_l","visits_past_year","medication_count","symptom_severity_l1_l4","functional_limitation_l1_l5","care_setting","smoking_status","diabetes","biomarker_positive"
"P001",62.9,25.6,122.1,2.47,8,4,"L4","L3","urgent","former","No","Yes"
"P002",48.9,25.6,120.5,1.25,1,,"L1","L2","outpatient","never","Yes","Yes"
```

Follow these rules when replacing the example with your data:

- Put column names in the first row and exactly one patient in every later row.
- Keep a complete, unique identifier such as `patient_id`. It is used only to
  label results and is excluded from clustering with `--id-column patient_id`.
- Put each candidate clustering feature in its own column. You may add, remove,
  or rename feature columns; they do not need to match the example names.
- Store continuous measurements as plain numbers without units or symbols;
  store count variables as non-negative integers.
- Store nominal categories as consistent text labels and binary variables as
  two consistent values such as `Yes` and `No`.
- Encode ordered categories as `L1`, `L2`, ... in increasing order when relying
  on automatic type inference. For other encodings, use a reviewed feature
  dictionary as described below.
- Represent missing values with an empty field, as in the second example row.
  `NA`, `N/A`, and `.` are also recognized as missing.
- Do not put a title, notes, units row, formulas, merged cells, or multiple tables
  in the CSV. Use UTF-8 text and comma separators.
- Supply at least 20 patients. A few hundred rows, like the intended use case,
  is preferable for stability assessment.
- Do not include cluster labels, outcomes, treatment arms, site, batch, or other
  variables that should not define clusters as ordinary features. Remove them or
  list them with `--exclude`, for example `--exclude outcome,site,batch`.

The example feature types are:

| Columns | Intended type |
|---|---|
| `age_years`, `bmi_kg_m2`, `systolic_bp_mmhg`, `crp_mg_l` | Continuous |
| `visits_past_year`, `medication_count` | Discrete counts |
| `symptom_severity_l1_l4`, `functional_limitation_l1_l5` | Ordinal (`L1` = lowest) |
| `care_setting`, `smoking_status` | Categorical |
| `diabetes`, `biomarker_positive` | Binary |

After the first run, inspect
`output/my_patients/feature_dictionary_used.csv`. If any inferred type or domain
is wrong, edit that dictionary and rerun:

```bash
Rscript analyze_csv.R \
  --input path/to/your_patients.csv \
  --dictionary output/my_patients/feature_dictionary_used.csv \
  --output-dir output/my_patients_reviewed \
  --pdf your_patients_reviewed_clustering_report.pdf \
  --id-column patient_id \
  --k auto
```

## Included mock CSV

The ready-to-use mock dataset is located at:

```text
data/mock_patients.csv
```

It contains 400 patient rows and 103 columns: 100 clustering features plus three
metadata columns. `patient_id` is the row identifier, `true_cluster` contains the
four simulation labels for evaluation only, and `site` is an excluded process
variable. Missing values are represented by empty CSV fields.

This simulation file is useful for benchmarking because its true clusters are
known. For the structure of an ordinary unlabeled input file, use
`data/example_patients_input.csv` instead.

Run the complete analysis directly against this CSV with:

```bash
Rscript analyze_csv.R \
  --input data/mock_patients.csv \
  --output-dir output/mock_csv_run \
  --pdf mixed_clustering_benchmark.pdf \
  --id-column patient_id \
  --truth-column true_cluster \
  --exclude site \
  --k auto
```

The command infers the mixed feature types from the CSV, runs every clustering
method, and generates the heatmap report. The truth and site columns are not
passed to imputation, dimension reduction, or clustering.

## Choosing the number of clusters

The default is `--k auto`; no cluster count is assumed. The workflow scans every
candidate from K=2 through `--max-k` (8 by default) and ranks candidates using
three unsupervised diagnostics: Gower/PAM silhouette, mixed-embedding k-means
silhouette, and latent-mixture BIC. Supplied truth labels are excluded from these
calculations and are used only for evaluation after clustering.

Adjust `--max-k` when a different search range is scientifically appropriate.
You can still provide an integer with `--k` for a planned sensitivity analysis,
but neither the ordinary CSV run nor `run_all.R` does so by default. The complete
ranking is saved as `auto_k_diagnostics.csv` in the chosen run directory.

To regenerate the deterministic mock CSV and its reference dictionary:

```bash
Rscript R/01_generate_mock_data.R
```

The reference feature dictionary is `data/feature_dictionary.csv`.

## Scripts

- `analyze_csv.R` is the one-command entry point for arbitrary CSV files.
- `R/csv_interface.R` validates the CSV, infers or reads the dictionary, selects
  K when requested, runs the benchmark, and renders the report.
- `R/01_generate_mock_data.R` generates the mock cohort, its complete-data audit
  copy, a feature dictionary, and a missingness summary.
- `R/clustering_functions.R` implements mixed-type preprocessing, stochastic
  hot-deck multiple imputation, all clustering algorithms, consensus clustering,
  label-invariant metrics, and MNAR perturbations.
- `R/02_run_clustering_benchmark.R` runs all methods over imputations and repeated
  80% subsamples, evaluates candidate values of K, performs MNAR sensitivity
  analysis, and saves CSV/RDS results.
- `R/03_render_report.R` renders the multipage PDF, including a dedicated dual-
  heatmap page for every method.
- `run_all.R` sends the generated mock CSV through the same public CSV interface.

## Methods implemented

1. Weighted Gower distance with partitioning around medoids (PAM)
2. Weighted Gower distance with average-linkage hierarchical clustering
3. k-prototypes
4. A dependency-light FAMD/PCAmix-style embedding followed by k-means
5. A diagonal Gaussian mixture on the mixed-data latent embedding
6. A typed latent class/profile mixture
7. Spectral clustering from a Gower affinity matrix
8. Consensus clustering across methods and imputations

The simulation's true labels are used only for evaluation. They are not supplied
to any clustering or imputation function.

## Using your own data

At minimum, place one patient per row and run:

```bash
Rscript analyze_csv.R \
  --input path/to/patients.csv \
  --output-dir output/my_run \
  --k auto \
  --id-column patient_id \
  --exclude outcome,site,batch
```

If known labels are available strictly for evaluation, add
`--truth-column known_group`. That column is excluded from fitting.

Automatic inference distinguishes continuous, low-cardinality integer,
categorical, binary, and common ordered-label columns. Because numeric ordinal
variables cannot be identified reliably from values alone, review
`feature_dictionary_used.csv` and rerun with:

```bash
Rscript analyze_csv.R \
  --input path/to/patients.csv \
  --dictionary output/my_run/feature_dictionary_used.csv \
  --output-dir output/my_reviewed_run \
  --k auto
```

The dictionary has these required columns:

| Column | Meaning |
|---|---|
| `feature` | Exact data-column name |
| `type` | `continuous`, `discrete`, `ordinal`, `categorical`, or `binary` |
| `domain` | Clinical domain used to balance total Gower weight |

Optional dictionary columns include `role` for reporting and `mnar_direction`
(`-1`, `0`, or `1`) for directional sensitivity analysis.

The same interface is available from an R session:

```r
source(file.path("R", "csv_interface.R"))

analyze_csv(
  input_path = "path/to/patients.csv",
  dictionary_path = NULL,
  output_dir = "output/my_cohort",
  id_column = "patient_id",
  exclude = c("outcome", "site"),
  seed = 20260902,
  k = "auto",
  m_imputations = 20,
  n_subsamples = 50
)
```

If truth labels are absent, ARI/NMI/accuracy are reported as missing; internal,
stability, uncertainty, heatmaps, and sensitivity analyses still run.

Important adaptations for a real study include clinically reviewing feature and
domain weights, reconstructing non-alphabetic ordinal levels explicitly, choosing
an imputation model compatible with the substantive analysis, and evaluating a
range of cluster counts rather than fixing K from the simulation.

## Machine-readable outputs

Each chosen output directory contains the dictionary used, automatic-K
diagnostics, a run manifest, and a `results/` directory containing:

- `method_metrics.csv`
- `cluster_assignments.csv`
- `membership_certainty.csv`
- `k_selection.csv`
- `mnar_sensitivity.csv`
- `missingness_summary.csv`
- `subsample_stability.csv`
- `analysis_results.rds`

PDF reports are always written separately to `output/pdf/`.

The mock data include a complete-data truth copy for simulation auditing. The
benchmark reads only the observed dataset and feature dictionary.
