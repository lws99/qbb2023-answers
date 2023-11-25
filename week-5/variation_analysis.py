#!/usr/bin/env python

import sys
import numpy
import matplotlib.pyplot as plt
from tqdm import tqdm


f=open(sys.argv[1])



read_depth=[]
genotype_qual=[]
allele_freq=[]
all_effects={}

#./variation_analysis.py annotated_decomp_filt_A01.vcf 


for line in f:
    if line.startswith('#'):
        continue
    fields = line.split()
    info=fields[7].split(";")
    for i in info:
        symbol,value=i.split("=")
        if symbol == "AF":
            af=[float(x) for x in list(value.split(","))] 
        elif symbol=="ANN":
            ann=value.split("|")
            effects=[]
            counter=1
            while counter < len(ann):
                effects.append(ann[counter])
                counter += 15  #moves to the next ann chunck of info
                #print(counter)
    #print(fields[9:])

    qual=[(x.split(":")[1]) for x in fields[9:]]
    genotype_qual += qual
    depth=[(x.split(":")[2]) for x in fields[9:]]
    read_depth += depth
    allele_freq += af
    for e in effects:
        if e not in all_effects:
            all_effects[e] = 0
        all_effects[e] += 1 
        #print(fields)



#=====================================================================================================================================================
# Step 3.1: Read depth distribution

# Plot a histogram showing the distribution of read depth at each variant across all samples (e.g. if you had 10 variants and 5 samples, you’d have 50 data points).

# This information can be found in the sample specific FORMAT fields and the end of each line. Check the file header to decide which ID is appropriate.

# Make sure you label the panel appropriately.


#==========================================================================================================================================================
#preparing the extracted data for graphing (clean the data and convert to integer)
read_depth_cleaned = [[num for num in sublist if num != '.'] for sublist in read_depth]
flat_read_depth_cleaned=[count for sublist in read_depth_cleaned for count in sublist]
#print(flat_read_depth_cleaned)
for i in flat_read_depth_cleaned:
    i=int(i)
    flat_read_depth_cleaned[i]=i


genotype_qual_cleaned = [[num for num in sublist if num != '.' and num != 'e' and num != "-"] for sublist in genotype_qual]
flat_genotype_qual_cleaned=[count for sublist in genotype_qual_cleaned for count in sublist]
#print(flat_genotype_qual_cleaned)
for i in flat_genotype_qual_cleaned:
    i=int(i)
    flat_genotype_qual_cleaned[i]=i

#===========================================================================================================================================================

#Graphing the data 

figure, ax = plt.subplots(ncols=2, nrows=2)
ax[0,0].set_xlabel("Sample")
ax[0,0].set_ylabel("Read Depth Frequency")
ax[0,0].set_title("Distribution of Read Depth at Each Variant")
ax[0,0].hist(flat_read_depth_cleaned) 





ax[0,1].set_xlabel("Sample")
ax[0,1].set_ylabel("Genotype Quality Frequency")
ax[0,1].set_title("Distribution of Genotyping Quality at Each Variant")
ax[0,1].hist(flat_genotype_qual_cleaned)





ax[1,0].set_xlabel("Allele Frequency")
ax[1,0].set_ylabel("Frequency")
ax[1,0].set_title("Allele Frequency of Variants")
ax[1,0].hist(allele_freq)




vals = sorted(all_effects.items(), key=lambda x:x[1], reverse=True)
labels, vals = zip(*vals[:10])
x = numpy.arange(len(labels))
ax[1,1].bar(x, vals)
ax[1,1].set_xlabel("Predicted Effects")
ax[1,1].set_ylabel("Frequency")
ax[1,1].set_title("Predicted Effects of Variants")




plt.show()


fig.savefig("data_exploration_week5.png")


plt.close(figure)




















# ['chrI', '968', '.', 'C', 'T', '40.1889', '.', 'AB=0.16129;ABP=33.9012;AC=8;AF=0.4;AN=20;AO=10;CIGAR=1X1M1X;DP=120;DPB=120;DPRA=0.428571;EPP=24.725;
# EPPR=3.0103;GTI=4;LEN=1;MEANALT=1;MQM=53.8;MQMR=52.0727;NS=10;NUMALT=1;ODDS=1.50578;PAIRED=0;
# PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=355;QR=3440;RO=110;RPL=0;RPP=24.725;RPPR=3.0103;RPR=10;RUN=1;SAF=10;SAP=24.725;SAR=0;SRF=110;SRP=241.872;SRR=0;TYPE=snp;
# ANN=T|upstream_gene_variant|MODIFIER|YAL067W-A|YAL067W-A|transcript|YAL067W-A_mRNA|protein_coding||c.-1512C>T|||||1512|,T|downstream_gene_variant|MODIFIER
# |YAL069W|YAL069W|transcript|YAL069W_mRNA|protein_coding||c.*319C>T|||||319|,T|downstream_gene_variant|MODIFIER|YAL068W-A|YAL068W-A|transcript|YAL068W-A_mRNA
# |protein_coding||c.*176C>T|||||176|,T|downstream_gene_variant|MODIFIER|PAU8|YAL068C|transcript|YAL068C_mRNA|protein_coding||c.*839G>A|||||839|,T|intergenic_region
# |MODIFIER|YAL068W-A-PAU8|YAL068W-A-YAL068C|intergenic_region|YAL068W-A-YAL068C|||n.968C>T||||||', 'GT:GQ:DP:AD:RO:QR:AO:QA:GL', 
# '0|1:0.00221008:18:16,2:16:482:2:67:-0.931828,0,-38.2477', '0|0:108.366:25:25,0:25:748:0:0:0,-7.52575,-62.6149', 
# '0|0:81.2777:16:16,0:16:459:0:0:0,-4.81648,-38.2989', '1|1:8.64119e-05:1:0,1:0:0:1:33:-3.29863,-0.30103,0', 
# '0|1:22.8698:13:10,3:10:372:3:111:-6.42438,0,-29.9006', '1|1:0.000463117:2:0,2:0:0:2:77:-7.28565,-0.60206,0', 
# '0|0:114.372:27:27,0:27:862:0:0:0,-8.12781,-73.3301', '0|0:51.1747:6:6,0:6:165:0:0:0,-1.80618,-7.32368', '1|1:0.000463117:2:0,2:0:0:2:67:-6.34564,-0.60206,0', 
# '0|0:63.2159:10:10,0:10:352:0:0:0,-3.0103,-31.9591']



#TRASH

# grab what you need from `fields`
#flattened_list = [count for sublist in read_depth_cleaned for count in sublist] #add a title to the plot numpy log10 conversion