suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(shinycssloaders)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(mongolite)
  library(glue)
  library(ggplot2)
  library(parallel)
  library(yaml)
})

freq_distribution <- function(n, lower_bound=0.0001, mean=0.10, sd=0.045) {
  pmax(lower_bound, rnorm(n, mean=mean, sd=sd))
}

datapath <- "pbmc_2018-05-24_dev"

genes_df <- readRDS(file.path(datapath, "genes.rds"))
i_df     <- readRDS(file.path(datapath, "individuals.rds"))

build_for_individual <- function(i_id, genes_df, output_dir) {
  idf <- genes_df %>% mutate(individual_id = i_id,
                             cat5 = sample(0:20, size=nrow(genes_df), replace=TRUE),
                             cat4 = sample(0:20, size=nrow(genes_df), replace=TRUE),
                             cat3 = sample(0:10, size=nrow(genes_df), replace=TRUE),
                             cat2 = sample(0:3, prob=c(5, 1, 1, 1), size=nrow(genes_df), replace=TRUE),
                             cat1 = sample(0:1, prob=c(5, 1), size=nrow(genes_df), replace=TRUE),
                             cat5maf = freq_distribution(nrow(genes_df)),
                             cat4maf = freq_distribution(nrow(genes_df)),
                             cat3maf = freq_distribution(nrow(genes_df)),
                             cat2maf = freq_distribution(nrow(genes_df)),
                             cat1maf = freq_distribution(nrow(genes_df)))
  message("Saving: ", i_id)
  saveRDS(idf, file=file.path(datapath, "variants_by_gene", glue("individual_{i_id}.rds")))
  TRUE
}

nothing <- i_df$individual_id %>%
  lapply(build_for_individual, genes_df)

