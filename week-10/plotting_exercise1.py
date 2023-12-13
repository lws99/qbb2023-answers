#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# normalize
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

# log
counts_df_logged = np.log2(counts_df_normed + 1)

# merge with metadata
full_design_df = pd.concat([counts_df_logged, metadata], axis=1)



#===================================================================================================================================================================
# Step 1.1: Distribution of expression across genes

# For subject GTEX-113JC, plot the distribution of expression (logged normalized counts) across all genes, excluding any genes with 0 counts. 
# Upload this figure for the assignment.


GTEX113JC_exp= counts_df_logged.loc["GTEX-113JC"]


GTEX113JC_exp_clean = GTEX113JC_exp[GTEX113JC_exp != 0]
#print(GTEX113JC_exp_clean)

fig, ax = plt.subplots()
ax.hist(GTEX113JC_exp_clean)
ax.set_title("Distribution of logged normalized counts, Subject GTEX-113JC")
ax.set_xlabel("Genes")
ax.set_ylabel("Counts")

#plt.show()


#===============================================================================================================================================
#1= male 2=female
    
# Step 1.2: Expression of a single gene between sexes

# For the gene MXD4, plot the distribution of expression (logged normalized counts) in males versus females. 
# Upload this figure for the assignment.

concat_expr_MD=pd.concat([metadata, counts_df_logged], axis=1, join="inner")
#print(concat_expr_MD.loc[:,"MXD4"])

gene_MXD4_expr=concat_expr_MD.loc[:,"MXD4"]
gene_MXD4_sexes=concat_expr_MD.loc[:,"SEX"]
#print(gene_MDX4_sexes)

MXD4_sex_expr=pd.concat([gene_MXD4_sexes, gene_MXD4_expr], axis=1, join="inner")
#print(MXD4_sex_expr)

MXD4_F_expr=MXD4_sex_expr[MXD4_sex_expr["SEX"]==2]
#print(MXD4_F_expr)

MXD4_M_expr=MXD4_sex_expr[MXD4_sex_expr["SEX"]==1]
#print(MXD4_M_expr)


#print(MXD4_M_expr)

fig2, ax = plt.subplots()

ax.hist(MXD4_M_expr["MXD4"], color = 'blue', alpha=0.7, label="Male")
ax.hist(MXD4_F_expr["MXD4"], color ='red', alpha=0.7, label="Female")
ax.set_title("Distribution of logged normalized counts for MXD4")
ax.set_xlabel("Subject") #not sure what this should be 
ax.set_ylabel("Counts")
plt.legend()



#===============================================================================================================================================================

# Step 1.3: Distribution of subject ages

# Plot the number of subjects in each age category. Upload this figure for the assignment.


age=concat_expr_MD.loc[:,"AGE"]
print(age)

fig3, ax = plt.subplots()


#FIX THE ORDER OF THE AGES TO BE LOWEST TO HIGHEST
ax.hist(age, bins=6)
ax.set_title("Number of subjects in each age group")
ax.set_xlabel("Age")
ax.set_ylabel("Count of Subjects")


#================================================================================================================================================================

# Step 1.4: Sex-stratified expression with age

# For the gene LPXN, plot the median expression (logged normalized counts) over time (i.e. in each age category), stratified by sex. Upload this figure for the assignment.


#stratified means make one line plot for females one for males, making a line plot, time on x axis


LPXN_expr=concat_expr_MD.loc[:,"LPXN"]
#print(gene_LPXN)
LPXN_sex=concat_expr_MD.loc[:,"SEX"]
LPXN_age=concat_expr_MD.loc[:,"AGE"]
LPXN_all_data=pd.concat([LPXN_age, LPXN_sex, LPXN_expr], axis=1, join="inner")
#print(LPXN_all_data)

#separate by sex
LPXN_all_data_M=LPXN_all_data[LPXN_all_data["SEX"]==1]
LPXN_all_data_F=LPXN_all_data[LPXN_all_data["SEX"]==2]

#separate by age 
LPXN_M_age_60=LPXN_all_data_M[LPXN_all_data_M["AGE"]== "60-69"]
LPXN_F_age_60=LPXN_all_data_F[LPXN_all_data_F["AGE"]== "60-69"]

LPXN_M_age_50=LPXN_all_data_M[LPXN_all_data_M["AGE"]=="50-59"]
LPXN_F_age_50=LPXN_all_data_F[LPXN_all_data_F["AGE"]=="50-59"]

LPXN_M_age_30=LPXN_all_data_M[LPXN_all_data_M["AGE"]=="30-39"]
LPXN_F_age_30=LPXN_all_data_F[LPXN_all_data_F["AGE"]=="30-39"]

LPXN_M_age_20=LPXN_all_data_M[LPXN_all_data_M["AGE"]=="20-29"]
LPXN_F_age_20=LPXN_all_data_F[LPXN_all_data_F["AGE"]=="20-29"]
#print(LPXN_F_age_20)


#find the median expression of LPXN
median_LPXN_M_age_60=np.median(LPXN_M_age_60.loc[:,"LPXN"])
median_LPXN_F_age_60=np.median(LPXN_F_age_60.loc[:,"LPXN"])

median_LPXN_M_age_50=np.median(LPXN_M_age_50.loc[:,"LPXN"])
median_LPXN_F_age_50=np.median(LPXN_F_age_50.loc[:,"LPXN"])

median_LPXN_M_age_30=np.median(LPXN_M_age_30.loc[:,"LPXN"])
median_LPXN_F_age_30=np.median(LPXN_F_age_30.loc[:,"LPXN"])

median_LPXN_M_age_20=np.median(LPXN_M_age_20.loc[:,"LPXN"])
median_LPXN_F_age_20=np.median(LPXN_F_age_20.loc[:,"LPXN"])


maleMedians = [median_LPXN_M_age_20, median_LPXN_M_age_30, median_LPXN_M_age_50, median_LPXN_M_age_60]
female_medians=[median_LPXN_F_age_20, median_LPXN_F_age_30, median_LPXN_F_age_50, median_LPXN_F_age_60]

#plot the median expression


fig4, ax=plt.subplots()
labels=["20-29", "30-39", "50-59", "60-69"]
#x axis is age 
#y is median expression
ax.scatter(range(len(maleMedians)), maleMedians, label="Male")
ax.scatter(range(len(female_medians)), female_medians, label="Female")
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(labels=labels, rotation=90)
ax.legend()
#ax.barplot(median_LPXN_F_age_60)



#plt.show()



#in step 2 you need to download the csv file for the data tro the left of the info page 
















#choose 3 key points/patterns observed and make 3 diff plots with that data