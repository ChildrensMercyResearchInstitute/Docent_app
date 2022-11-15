#Docent app

## Data Prepration


###(a) Generate fake affected gene data
    Rscript generate_fake_affected_gene_data.r
###(b) Export Seruat object to Rds file.
    Rscript export_seurat_to_rds.r
###(c) Generate genotype.csv file from multisample vcf
    Rscript import_vcf_to_genotype_csv.r
###(d) Import genotype csv file into sqlite database
    Rscript import_genotypes_to_sqlite.sh


##Build and push docker image
sh docker_build.sh

##Running the app
docker run -p 3838:3838 -v /path to data_folder/pbmc_2018-05-24_dev:/docent_data docker.io/man4ish/docent:latest

path to data folder on gpfs : /home/mkumar1/app_backup/backup_20220708/Docent_app/pbmc_2018-05-24_dev

### Example
docker run -p 3838:3838 -v /home/mkumar1/app_backup/backup_20220708/Docent_app/pbmc_2018-05-24_dev:/docent_data docker.io/man4ish/docent:latest


##browse the application
http://0.0.0.0:3838
