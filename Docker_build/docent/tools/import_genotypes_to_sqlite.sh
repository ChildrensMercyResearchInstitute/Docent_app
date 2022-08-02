sqlite3 genotype.db
CREATE TABLE genotypes(i TEXT, rsid TEXT, geno TEXT);
CREATE INDEX geno_index ON genotypes (rsid);
.mode csv
.import genotypes.csv genotypes
