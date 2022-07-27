# Docent: a single-cell gene expression exploration tool

Docent is a web-based single-cell gene expression exploration tool that integrates individual genotype and phenotype data. 

## Requirements

 - [R](https://www.r-project.org) 3.5.0 or higher (TODO: list of required R packages)
 - A [Seurat](https://satijalab.org/seurat) object that has been saved to disk via `saveRDS()`
 - SQLite (https://www.sqlite.org/)
 - A VCF file providing genotypes for each individual

## Installation: Shiny Server

This tool is a Shiny application that can be hosted with Shiny Server by following these steps:

  1. Checkout the source repository
  2. Create a directory within the ShinyApps directory of the user that will be hosting the app (default `/home/username/ShinyApps`)
  3. Copy `app.R` and `config.yaml.sample` to the directory from step 2
  4. Rename `config.yaml.sample` to `config.yaml` and modify settings as needed (see Configuration below)
  5. Navigate to the Shiny Server's URL with the appropriate username and the directory name from step 2 (typically something like http://shiny-server:3838/username/app_directory)
  
## Configuration

The `config.yaml` file will be read by the app on startup. Copy the `config.yaml.example` to `config.yaml` and modify as needed:

  - `datapath`: Path to a data directory generated using the `export_seurat_to_rds.r` script
  - `individual_remapping`: Path to a saved R data.frame specifying any necessary mappings between individuals in the single-cell data and those in the genotype data
  - `cores`: Maximum number of threads to use during a single request
