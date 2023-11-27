Step 1.1 commands:

plink --noweb --vcf genotypes.vcf --pca 10




Step 2.1:

plink --noweb  --vcf genotypes.vcf --freq --out allele_frequencies



Step 3.1:

plink --noweb  --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex --out CB1908_IC50_GWAS

plink --noweb  --vcf genotypes.vcf --linear --pheno GS451_IC50.txt.txt --covar plink.eigenvec --allow-no-sex --out GS451_IC50_GWAS



Question 3.4:

