#!/usr/bin/env python




import subprocess
import matplotlib.pyplot as plt 
import numpy as np
#===========================================================================================================================================
# Exercise 1: Performing PCA
# Step 1.1: Compute genotype PCs

# First, you want to perform PCA to compute the top genotype PCs for the samples in your data set. This is useful not only for visualizing the population structure in your data set, but also allows you to control for this structure when you perform your GWAS later.

# Using plink, perform PCA on the genotypes of the cell lines in the data set. Output the top 10 principal components.

# Record the plink command(s) you used in your README.md.





# plink --noweb --vcf genotypes.vcf --pca 10

# head plink.eigenval principal components
# 2.23286
# 1.91716
# 1.88858
# 1.83419
# 1.82363
# 1.81699
# 1.80463
# 1.80033
# 1.7953
# 1.78938





PC=np.loadtxt("plink.eigenvec")
#print(plink)


PC1=PC[:, 2]

PC2=PC[:, 3]

fig, ax = plt.subplots()
ax.scatter(PC1, PC2, label="PC1 vs PC2") 
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PC1 vs PC2")
#plt.show()



#====================================================================================================================================================================
# Exercise 2: The allele frequency spectrum
# Step 2.1: Compute allele frequencies

# The shape of the allele frequency specturm (AFS) can tell you interesting things about your sample, such as demographic history, and evidence of selection (read more here. While you will not be running these specific tests today, it is still interesting (and good practice) to visualize the AFS.

# Usually, the AF of each variant is stored in the INFO field in your VCF. The vcf file we’re using does not have an INFO field, so you’ll need to calculate AF yourself.

# Using plink, calculate allele frequencies for the variants in your VCF and output these to a file.

# Record the plink command(s) you used in your README.md.


#plink --noweb  --vcf genotypes.vcf --freq --out allele_frequencies





#====================================================================================================================================================================
# Step 2.2: Plot AFS

# Now that you’ve computed allele frequencies, you can visualize the AFS of the samples in your data set.

# In your plotting.py script, plot the AFS of the samples in your data set as a histogram. Ensure that your plot is nicely formatted and properly labeled.

AF=np.loadtxt("allele_frequencies.frq", skiprows=1, usecols=4)
#print(AF)


fig, ax1 = plt.subplots()
ax1.hist(AF, bins=80) 
ax1.set_xlabel("Allele Frequency")
ax1.set_ylabel("Counts")
ax1.set_title("Distribution of Allele Frequencies")
#plt.show()


#======================================================================================================================================================================


# Step 3.1: Running the GWAS

# Now you’re ready to run your GWAS!

# Using plink, perform quantitative association testing for each of your two phenotypes (you will be running plink one each per phenotype). 
# Use the top 10 principal components you computed in Step 1.1 as covariates in your analysis to control for non-independence due to relatedness/population structure. 
# <li> Be sure to use the --allow-no-sex option</li> <li> You may find this portion of the plink documentation helpful for performing association testing on 
# each of the phenotypes.</li> <li> HINT: We do NOT have case/control phenotpe data, so --assoc is not the correct plink flag




#plink --noweb  --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex --out CB1908_IC50_GWAS

#plink --noweb  --vcf genotypes.vcf --linear --pheno GS451_IC50.txt.txt --covar plink.eigenvec --allow-no-sex --out GS451_IC50_GWAS




#=================================================================================================================================================================


# Step 3.2: Visualing GWAS results

# Now that you’ve run your GWAS analyses, you want to visualize the results.

# In your plotting.py script, produce a two-panel figure (one column, two rows) depicting the Manhattan plot for each of your two GWAS analyses 
# (i.e. each phenotype). Each Manhattan plot should be one of the two panels in your figure. In each Manhattan plot, highlight SNPs with p-values 
# less than 10-5 in a different color.

# Ensure that your figure/panels are nicely formatted and properly labeled.

#----------------------------------------------------------------------------------------

#load and clean the CB1908 data
p_val_list=[]
snp_pos_list=[]
effect_size_list=[]
for i,line in enumerate(open("CB1908_IC50_GWAS.assoc.linear")):
    if i == 0:
        continue 
    fields=line.strip().split()
    p_val=fields[8]
    p_val = float(p_val)
    p_val_list.append(p_val)
    chrom=int(fields[0])
    snp_pos=int(fields[2])
    snp_pos_list.append(snp_pos)
    snp_name=fields[1]
    effect_size=fields[6]
    effect_size_list.append(effect_size)


# # lowest pval is 8.199e-12
# #effect size for lowest pval is 2.002
   


#find -log10 of pvals of CB1908 data
log_p_val_list=[]
for i in p_val_list:
    log_p_val=-1*np.log10(i)
    log_p_val_list.append(log_p_val)

#find lowest pval of CB1908 data
min_CB_p_val_np=np.min(p_val_list)
#print(min_CB_p_val_np)



#organizing snp data by chromosome position - CB1908
pval_by_position={}
total_snp_count=0
for i,line in enumerate(open("CB1908_IC50_GWAS.assoc.linear")):
    if i == 0:
        continue 
    fields=line.strip().split()
    p_val=fields[8]
    p_val = float(p_val)
    snp_name=fields[1]
    if p_val==min_CB_p_val_np:
        top_snp_name=snp_name
    pval_by_position.setdefault(chrom, []) #don't put in a fake value if the chromosome doesn't exist
    pval_by_position[chrom].append((snp_pos, p_val)) #defined the values in the dictionary for each chrom key 
    total_snp_count += 1


#assign chromosome colors even chromosomes are blue odd are red - CB1908
colors=[]
pvalues=[]
x_position=[x for x in range(total_snp_count)]
for chrom in pval_by_position:
    color = ["blue", "red"][chrom % 2]
    sorted_pval=sorted(pval_by_position[chrom], key=lambda x: x[0]) #sort based on x where x is chromosome
    for i in sorted_pval:

        pvalues.append(i[1]) #pval info
        colors.append(color) #color info --> determined by what chromosome 



# # 12  rs10876043   49190411    G        ADD      161      2.002        7.422    8.199e-12


#-------------------------------------------------------------------------

#load in and clean GS451 data
p_val_list_GS=[]
snp_pos_list_GS=[]
effect_size_GS_list=[]
for i,line in enumerate(open("GS451_IC50_GWAS.assoc.linear")):
    if i == 0:
        continue 
    fields_GS=line.strip().split()
    p_val_GS=fields_GS[8]
    p_val_GS = float(p_val_GS)
    p_val_list_GS.append(p_val_GS)



    chrom_GS=int(fields_GS[0])

    snp_pos_GS=int(fields_GS[2])
    snp_pos_list_GS.append(snp_pos_GS)
    snp_name_GS=fields_GS[1]

    effect_size_GS=fields_GS[6]
    effect_size_GS_list.append(effect_size_GS)



   
#find -log10 of GS451 data
log_p_val_list_GS=[]
for i in p_val_list_GS:
    log_p_val_GS=-1*np.log10(i)
    log_p_val_list_GS.append(log_p_val_GS)

#find lowest pval for GS451 data
min_GS_p_val_np=np.min(p_val_list_GS)
#print(min_GS_p_val_np)




# small_p=[]
# for i in log_p_val_list_GS:
#     if i < 10e-5:
#         small_p.append(i)
# #print(small_p)



#create Manhattan plots for both datasets 

fig, ax2 = plt.subplots(nrows=2)
ax2[0].scatter(snp_pos_list, log_p_val_list)
ax2[0].set_title("Manhattan Plot for CB1908_IC50_GWAS")
ax2[0].set_xlabel("Chromosome Position")
ax2[0].set_ylabel("-log10 pval")

ax2[1].scatter(snp_pos_list_GS, log_p_val_list_GS)
ax2[1].set_title("Manhattan Plot for GS451_IC50_GWAS")
ax2[1].set_xlabel("Chromosome Position")
ax2[1].set_ylabel("-log10 pval")





#===============================================================================================================================================================

# # Step 3.3: Visualizing effect-size

# # You want to dig deeper into the top GWAS hits from your analyses.

# # Choose one of the traits for which you performed GWAS.

# # In your plotting.py script: for the top associated SNP (lowest p-value) of that trait, plot the effect size of that 
# # variant on the chosen trait by creating a boxplot of the phenotype stratified by genotype.


#matching genotypes and phenotypes

#find in genotypes.vcf for each individual their genotype, 3 lists of phenotypes sep by genotype 




for i,line in enumerate(open("genotypes.vcf")):
    if i<27:
        continue
    genotypes_strip_split=line.strip().split()
    if genotypes_strip_split[2] == top_snp_name:
        genotypes=genotypes_strip_split[9:]
    if i == 27:
        id_list=genotypes_strip_split[9:]
# print(genotypes)
# print(id_list)


heterozygous=[]
homozygous_1=[]
homozygous_0=[]

for i,line in enumerate(open("CB1908_IC50.txt")):
    if i == 0:
        continue
    phenotypes_strip_split=line.strip().split()
    phenotype_id= phenotypes_strip_split[0]+"_"+phenotypes_strip_split[1]
    #print(phenotype_id)
    for j,position in enumerate(id_list):
        if phenotype_id==position and phenotypes_strip_split[2] != "NA":
            if genotypes[j]=="0/0":
                homozygous_0.append(float(phenotypes_strip_split[2]))
            elif genotypes[j]=="0/1":
                heterozygous.append(float(phenotypes_strip_split[2]))
            elif genotypes[j]=="1/1":
                homozygous_1.append(float(phenotypes_strip_split[2]))
#print(homozygous_1, homozygous_0, heterozygous)





#plot these now 
fig, ax3 = plt.subplots()
ax3.boxplot([heterozygous, homozygous_1, homozygous_0])
ax3.set_title("Boxplot for CB1908_IC50_GWAS snp:rs10876043")
ax3.set_xlabel("Genotype")
ax3.set_ylabel("Effect Size of Variant")


plt.show()



















# #enumerate allows you to go line by line i is number of the line and line is the data in the line
















