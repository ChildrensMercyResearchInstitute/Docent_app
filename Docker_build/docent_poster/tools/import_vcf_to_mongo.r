library(VariantAnnotation)
library(dplyr)
library(magrittr)
library(readr)

vcf <- "/analysis/projects/xenon/impute2/vcf/All_Genotyping.vcf.gz"
x <- readGT(vcf)

sample_names <- colnames(x)
sample_names <- sample_names[!grepl("^NA.*", sample_names)]
names(sample_names) <- sample_names

genos_df <- sample_names %>%
  lapply(function(sample_name) {
    data_frame(rsid = rownames(x),
               genotype = x[, sample_name, drop=TRUE])
  }) %>%
  bind_rows(.id="individual") %>%
  filter(genotype != ".") %>%
  select(`i.string()`    = individual,
         `rsid.string()` = rsid,
         `geno.string()` = genotype)

write_csv(genos_df, "genotypes.csv")
