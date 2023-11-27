#!/usr/bin/env python


import sys
import numpy
import matplotlib.pyplot as plt

# Part 3a

# Using the above bedgraph files, write a Python script to perform a number of comparisons between the two sets of 
# methylation calls. Your script should be able to:

#comparing me calls between ONT (nanopore seq) and bisulfite (bismark)

# Parse the bedgraph files
# Calculate the number of sites present only in the bismark file, present only in the nanopore file, and the shared sites 
# as a percentage of total sites (both unique and shared sites) and record them in your README.md file.

#cpg file has one cpg location per line (how much me they found, one me= one line)

#./week-7_assign.py ONT.cpg.chr2.bedgraph normal.ONT.chr2.bedgraph tumor.ONT.chr2.bedgraph bisulfite.cpg.chr2.bedgraph normal.bisulfite.chr2.bedgraph tumor.bisulfite.chr2.bedgraph x
def main():
    # Load file names from command 
    cpg_ONT, ONT_normal, ONT_tumor, cpg_bisulf, normal_bisulf, tumor_bisulf, out_fname = sys.argv[1:8]

    # Load data from files
    bisulfite=load_data(cpg_bisulf)
    bisulfite_n= load_data(normal_bisulf)
    bisulfite_t=load_data(tumor_bisulf)
    ONT=load_data(cpg_ONT)
    ONT_n=load_data(ONT_normal)
    ONT_t=load_data(ONT_tumor)

    # Find reads that appear more than once in datasets (read appears more than once in output of mapping) --> set shows you where you have unique vs overlapping names
    bisulfite_set=set(bisulfite.keys())   #keys are from the load_data function --> cols 1 and 2 from the data (chr and starting position)
    ONT_set=set(ONT.keys())
    #print(type(ONT_set))
    #intersect takes only items that are in both dictionaries 
    ONT_bs_overlap=bisulfite_set.intersection(ONT_set)
    #union takes items in either or both dictionaries
    ONT_bs_together=bisulfite_set.union(ONT_set)


    #print(ONT_bs_overlap)

    #how many are only bisulfite reads
    print(f"bisulfite only reads: {len(bisulfite_set)-len(ONT_bs_overlap)}")
    BS_only=len(bisulfite_set)-len(ONT_bs_overlap)
    print(BS_only)

    #how many are ONT reads
    print(f"ONT only reads: {len(ONT_set)-len(ONT_bs_overlap)}")
    ONT_only=len(ONT_set)-len(ONT_bs_overlap)



    #how many reads are shared between ONT and bisulfite reads
    print(f"percentage of shared reads: {len(ONT_bs_overlap)/len(ONT_bs_together)*100} %")

    list_ONT= list(ONT_set)
    #print(list_ONT[:10])
    flattened_list_ONT=[]
    for i in list_ONT:
        for j in i:
            if j != "chr2":
                j=int(j)
                #print(type(j))
                flattened_list_ONT.append(j)
    #print(flattened_list_ONT)
    list_BS= list(bisulfite_set)
    #print(list_BS[:10])
    flattened_list_BS=[]
    for i in list_BS:
        for j in i:
            if j != "chr2":
                j=int(j)
                #print(type(j))
                flattened_list_BS.append(j)
    #print(flattened_list_BS[:10])




    #print(flattened_list_ONT[:10])

    fig, ax = plt.subplots()

    #I don't thiink these numbers are quite right --> need to be coverage and not the position of chromosome(the ONT only and BS only reads)
    ax.hist(flattened_list_ONT, color="red", alpha=0.7, label="nanopore")
    ax.hist(flattened_list_BS, color="blue", alpha=0.6, label="bisulfite")
    ax.set_title("")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Frequency")
    plt.show()
    





    # for i in range(len(normal)):         
    #     if normal[i][0] not in normal_set:               
    #         normal_set.add(normal[i][0])
    #     else:
    #         normal_multi.add(normal[i][0])
    # normal_single=normal_set.difference(normal_multi)


def load_data(fname):
    data={}
    for line in open(fname):
        line=line.rstrip().split("\t")   #\t split on tab
        key=(line[0], line[1])
        data[key]=(float(line[3]), int(line[4]))
    return(data)  







main()










# Part 3b:

# Plot the distribution of coverages across CpG sites for each track on the same plot. Make sure to indicate which distribution corresponds to which track. 
# In order to visualize both distributions on the same plot, it may be useful to use the alpha option to set the transparency of the plotted data.

# Q2: How does using nanopore for methylation calling differ from bisulfite sequencing in terms of coverage? Which method appears better and why?









