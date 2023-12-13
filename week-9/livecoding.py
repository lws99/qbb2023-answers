#!/usr/bin/env python

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest
from pydeseq2 import preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
from matplotlib import pyplot as plt




#whole blood (755 total individuals) --> looking at every gene expressed in blood of males and females in a population


#in gtex_metadata.txt :
#1 is males 2 is females (I think)
#DTHHRDY--> cause of death on the hardy scale (ex: 0 is on a ventilator)


#=====================================================================================================================================================================

# Step 1.1: Loading data and importing libraries

# Because the data are rectangular data of mixed type, you can use pandas to read the data into a data frame. You will also need statsmodels to 
# perform your “homemade” differential expression test, followed by PyDESeq2 to apply the more sophisticated test. You’ll also need numpy to do some data transformations. 
# Here is code to load those up, along with the data. You will need to modify the file paths to direct it to wherever you stored the data.


# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

#print(counts_df)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)




#correct count matrix for our differential analysis matrix --> to give you a normalized matrix (normalizing for the gene expression and )


counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

counts_df_normed = np.log2(counts_df_normed + 1)

gene_list=counts_df_normed.columns.tolist()
gene_count=len(gene_list)
#print(gene_list[1:10])

#bind together (concatinate) the different matricies together (make a big dataframe)


full_design_df = pd.concat([counts_df_normed, metadata], axis=1)
#print(full_design_df)



#fitting linear models going gene by gene (looping over column 2 of the full_design_df --> looping over each gene and determining if there are differences in expression
#based on sex)
#intercept= average expression overall --> we don't care abut this WE CARE ABOUT IF THE SLOPE IS SIGNIFICANTLY DIFFERENT FROM 0 

model = smf.ols(formula = 'Q("DDX11L1") ~ SEX', data=full_design_df)
results = model.fit()
#print(results.summary())


#This gene is not differentially expressed 
#SEX           expression=-0.0816      0.076     -1.075      p=0.283      -0.231       0.068 expression is not siginficantly different from 0 based on p 


#grabbing the slope and the pvalue 
slope = results.params[1]
pval = results.pvalues[1]




#====================================================================================================================================================================


# Step 1.5: Extend this test to all genes

# Write a for-loop in Python to extend this test to all genes in your matrix. For each gene that you test, store the slopes and p-values, along 
# with the associated gene names in a useful data structure (pandas DataFrame, arrays, …, ?).

# NOTE: This step will probably take a few minutes to run. We recommend outputing your results to a file, and then reading it back in to your 
# script so that you don’t have to re-run this whole process each time you run your python script.

# After running the analysis on all the genes, use the Benjamini-Hochberg procedure to determine which genes are sex-differentially expressed at an FDR of 10%. 
# Statsmodels has a great function to do this correction for you: statsmodels.stats.multitest.fdrcorrection

# Write the list of genes that are differentially expressed at a 10% FDR to a file. You will be uploading this file with your assignment.

#genes are the columns



#output_writing=open("output_file.txt", "w")
# print(genes)

path = "output_file_all_data.txt"

if not os.path.exists(path):
    all_genes_test=pd.DataFrame({"slope":pd.Series(float), "pval":pd.Series(float)}, index=gene_list)

    for gene in gene_list:
        #gene=full_design_df.loc[:,1:]
        all_model = smf.ols(formula = f'Q("{gene}") ~ SEX', data=full_design_df)
        all_results = all_model.fit() 
        all_slope = all_results.params[1]
        all_genes_test.loc[gene, "slope"] = all_slope
        all_pval = all_results.pvalues[1]
        all_genes_test.loc[gene, "pval"] = all_pval
    
    #get adjusted pvals at 10% FDR
    all_genes_test["pval"].fillna(1) #all NA's will now be one, will not include pval becuase it is 1 
    all_genes_test["adjusted_pval"]=multitest.fdrcorrection(all_genes_test["pval"], alpha=0.1, method='indep')[1] #indep bc expecting to be independent false[1] only takes the adjusted pval from the fdr
    all_genes_test["significant_pval"]=all_genes_test["adjusted_pval"] < 0.1
    all_genes_test.to_csv("output_file_all_data.txt", header=True, index=True, sep='\t')  #when I open it again slope and pval need to be floats not strings


#open the file
all_genes_output_file=pd.read_csv("output_file_all_data.txt", header=0, index_col=0, sep='\t')

#grab diff expressed genes
diff_expr_genes=all_genes_output_file.loc[all_genes_output_file["significant_pval"]].index.tolist()

#write diff expr genes to a file
with open("diff_expr_genes.txt", "w") as file:
    for gene in diff_expr_genes:
        file.write(f"{gene}\n")

#exercise 2============================================================


dseq_path="dseq_results.txt"
if not os.path.exists(dseq_path):
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="SEX",
        n_cpus=4,
    )


    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()
    results_dds = stat_res.results_df
    results_dds["padj"].fillna(1)
    results_dds["significant_pval"]=results_dds["padj"] < 0.1
    #write dseq to file
    results_dds.to_csv("dseq_results.txt", header=True, index=True, sep='\t') 

#open the file
dseq_output_file=pd.read_csv("dseq_results.txt", header=0, index_col=0, sep='\t')

#get diff expr genes dseq
diff_expr_genes_dseq=dseq_output_file.loc[dseq_output_file["significant_pval"]].index.tolist()

#write diff expr genes to a file
with open("d_seq_diff_expr_genes.txt", "w") as file:
    for gene in diff_expr_genes_dseq:
        file.write(f"{gene}\n")





# Compare the list of genes that ARE differentially expressed at a 10% FDR to those you identified in Step 1.5. What is the percentage of overlap? 
# Compute this percentage in your code as a “Jaccard index”, which is defined as the intersection divided by the union: 
# ((number of genes that were significant in steps 1 and 2) / (number of genes that were significant in steps 1 or 2)) * 100%. 
# Record this in your README.md file for this assignment, which you will upload with the assignment.

# gene_intersection=diff_expr_genes_dseq.intersection(diff_expr_genes)
# gene_union=diff_expr_genes_dseq.union(diff_expr_genes)



# def (significant_genes, significant_genes_dseq):
#     Jaccard_index=(len(significant_genes[:,0])/len(significant_genes_dseq))*100
#     return(Jaccard_index)


def Jaccard(manual, dseq): #set1= manual, set2=dseq
    Jaccard_index=(len(set(manual)&set(dseq))/len(set(manual)|set(dseq))) *100
    return(Jaccard_index)

genes_Jaccard=Jaccard(diff_expr_genes, diff_expr_genes_dseq)
print(genes_Jaccard)



#Exercise 3 ======================================================================

fig, ax = plt.subplots()

log2_fold_change=dseq_output_file.loc[:,"log2FoldChange"]
pval=dseq_output_file.loc[:,"padj"]
log_pval=np.log10(pval) * -1

ax.scatter(log2_fold_change, log_pval)
ax.set_xlabel("log2FoldChange")
ax.set_ylabel("-log10 pval")
ax.set_title("PyDESeq2 Differential Expression Results")\



#use this to change the color of the datapoints
ax.scatter([0.2, 0.5], [75, 99])



plt.show()


    







