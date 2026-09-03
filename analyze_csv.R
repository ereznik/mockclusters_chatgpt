#!/usr/bin/env Rscript

# One-command entry point for any patient-by-feature CSV.
# Run from the project directory; use --help to see all options.
source(file.path("R", "csv_interface.R"))
run_csv_cli()
