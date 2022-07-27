sqlite3 genotype.db
CREATE TABLE genotypes(i TEXT, rsid TEXT, geno TEXT);
CREATE INDEX geno_index ON genotypes (rsid);
.mode csv
.import /home/mkumar1/ShinyApps/docent_poster/genotypes.csv genotypes
