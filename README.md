# Mixed-type patient clustering benchmark

This project accepts any patient-by-feature CSV, executes eight mixed-data
clustering strategies, and creates a heatmap-heavy PDF report. It also includes a
400-patient, 100-feature mock cohort with four known clusters, correlated noise,
outliers, and structured missingness.

## Quick start

Analyze any CSV from this directory:

```bash
Rscript analyze_csv.R --input path/to/patients.csv --k auto
```

The command infers a feature dictionary if none is supplied and writes it to the
run directory for review. Use `Rscript analyze_csv.R --help` for every option.

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

## Included mock CSV

The ready-to-use mock dataset is located at:

```text
data/mock_patients.csv
```

It contains 400 patient rows and 103 columns: 100 clustering features plus three
metadata columns. `patient_id` is the row identifier, `true_cluster` contains the
four simulation labels for evaluation only, and `site` is an excluded process
variable. Missing values are represented by empty CSV fields.

Run the complete analysis directly against this CSV with:

```bash
Rscript analyze_csv.R \
  --input data/mock_patients.csv \
  --output-dir output/mock_csv_run \
  --pdf output/pdf/mixed_clustering_benchmark.pdf \
  --id-column patient_id \
  --truth-column true_cluster \
  --exclude site \
  --k 4
```

The command infers the mixed feature types from the CSV, runs every clustering
method, and generates the heatmap report. The truth and site columns are not
passed to imputation, dimension reduction, or clustering.

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
diagnostics, a run manifest, the PDF (unless directed elsewhere), and a `results/`
directory containing:

- `method_metrics.csv`
- `cluster_assignments.csv`
- `membership_certainty.csv`
- `k_selection.csv`
- `mnar_sensitivity.csv`
- `missingness_summary.csv`
- `subsample_stability.csv`
- `analysis_results.rds`

The mock data include a complete-data truth copy for simulation auditing. The
benchmark reads only the observed dataset and feature dictionary.
