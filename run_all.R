# Reproduce the complete CSV-first simulation, benchmark, and PDF report.
source(file.path("R", "01_generate_mock_data.R"))
generate_mock_patients(n = 400L, seed = 20260902L, output_dir = "data")

source(file.path("R", "csv_interface.R"))
analyze_csv(
  input_path = file.path("data", "mock_patients.csv"),
  output_dir = file.path("output", "mock_csv_run"),
  pdf_path = file.path("output", "pdf", "mixed_clustering_benchmark.pdf"),
  dictionary_path = NULL,
  id_column = "patient_id",
  truth_column = "true_cluster",
  exclude = "site",
  k = 4L,
  max_k = 8L,
  m_imputations = 8L,
  n_subsamples = 10L,
  seed = 20260902L
)
