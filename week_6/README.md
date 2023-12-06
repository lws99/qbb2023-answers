Step 1.1 commands:

plink --noweb --vcf genotypes.vcf --pca 10




Step 2.1:

plink --noweb  --vcf genotypes.vcf --freq --out allele_frequencies



Step 3.1:

plink --noweb  --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex --out CB1908_IC50_GWAS

plink --noweb  --vcf genotypes.vcf --linear --pheno GS451_IC50.txt.txt --covar plink.eigenvec --allow-no-sex --out GS451_IC50_GWAS



Question 3.4:

DIP2B is the closest gene to the highest associated SNP on chromosome 12 for CB1908_IC5. This gene is required for normal synaptic transmission, so changes in this gene could possibly cause neurodegenerative diseases. In this case, this SNP could be affecting lymphocyte generation becuase the signaling for lymphocyte synthesis could be interrupted. 

SHC2 is the closest gene to the highest associated SNP on chromosome 19 for CB1908_IC50. This gene is required for coupling activated growth factor receptors to signaling pathways in neurons. Changes in this gene could possibly cause cancer. In this case, this SNP could be affecting lymphocyte generation becuase growth factor receptors are unable to communicate with the body, therefore if more lymphocytes need to be make, the growth factor may not be able to signal this. 



