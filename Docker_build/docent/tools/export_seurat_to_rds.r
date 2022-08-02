suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-s", "--seurat"), 
              type="character",
              default=NA,
              help="Path to Seurat object"),
  make_option(c("-o", "--output"),
              type="character",
              default=NA,
              help="Path to output folder"),
  make_option(c("-i", "--individual"),
              type="character",
              default="individual",
              help="Seurat metadata column specifying individual"),
  make_option(c("-c", "--cluster"),
              type="character",
              default="cluster",
              help="Seurat metadata column specifying cluster"),
  make_option(c("-t", "--threads"),
              type="integer",
              default=12,
              help="Number of threads to use (default: 12)")
)

opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(file.exists(opt$seurat))

if(!dir.exists(opt$output)) {
  dir.create(opt$output)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(glue)
  library(parallel)
})


####################################

pn <- function(i, ...) {
  prettyNum(i, big.mark=",", ...)
}

faster_load <- function(file) {
  stopifnot(file.exists(file))
  gzip_in <- pipe(paste0("gzip -dc ", file), "rb")
  on.exit(close(gzip_in))
  get(load(gzip_in))
}

faster_readRDS <- function(file) {
  stopifnot(file.exists(file))
  gzip_in <- pipe(paste0("gzip -dc ", file), "rb")
  on.exit(close(gzip_in))
  readRDS(gzip_in)
}

####################################

message("Loading Seurat object...")
so <- faster_load(opt$seurat)

metadata_columns <- names(so@meta.data)

if(!opt$individual %in% metadata_columns) {
  message("Available metadata columns: ", paste0(metadata_columns, collapse=", "))
  stop(paste0("Individual metadata column not found: ", opt$individual))
}


metadata <- so@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode, 
         cluster      = !!opt$cluster, # tsne.default.ident, 
         individual   = !!opt$individual, # geno.name
         mito_percent = percent.mito,
         umi_count    = nUMI) %>%
  mutate(barcode_id = 1:length(barcode))

exp_data <- so@data
rm(so)
invisible(gc())

build_gene_expression <- function(gene_id, expression_matrix, metadata, gene_id_table) {
  gene_name <- gene_id_table$gene[gene_id_table$gene_id == gene_id]
  cell_exp <- expression_matrix[gene_name, ]
  data_frame(barcode = names(cell_exp), 
             gene_id = gene_id,
             exp = cell_exp) %>%
    inner_join(select(metadata, barcode, barcode_id), by="barcode") %>%
    select(-barcode)
}

gene_table <- data_frame(gene = rownames(exp_data),
                         gene_id = 1:length(gene))

individual_table <- data_frame(individual = unique(metadata$individual),
                               individual_id = 1:length(individual))

metadata %<>% inner_join(individual_table, by="individual") %>%
  select(barcode, barcode_id, cluster, individual_id)

clustering_dir <- file.path(opt$output, "clustering")
if(!dir.exists(clustering_dir)) dir.create(clustering_dir, recursive=TRUE)

saveRDS(metadata, file.path(clustering_dir, "default.rds"))
saveRDS(gene_table, file.path(opt$output, "genes.rds"))
saveRDS(individual_table, file.path(opt$output, "individuals.rds"))

genes <- gene_table$gene_id

message(glue("Exporting {pn(length(genes))} genes"))
nothing <- genes %>%
  mclapply(function(single_gene_id) {
    gene_exp <- build_gene_expression(single_gene_id, exp_data, metadata, gene_table)
    message("Writing: ", single_gene_id)
    output_path <- file.path(opt$output, "expression", single_gene_id %% 100)
    if(!dir.exists(output_path)) dir.create(output_path, recursive=TRUE)
    saveRDS(gene_exp, file.path(output_path, glue("gene_{single_gene_id}.rds")))
    invisible(gc())
  }, mc.cores=opt$threads)

message("Output complete.")

