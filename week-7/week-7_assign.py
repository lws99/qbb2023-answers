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


#usage:
#./week-7_assign.py ONT.cpg.chr2.bedgraph normal.ONT.chr2.bedgraph tumor.ONT.chr2.bedgraph bisulfite.cpg.chr2.bedgraph normal.bisulfite.chr2.bedgraph tumor.bisulfite.chr2.bedgraph x
#write out what each of the args is as well


def main():
    # Load file names from command 
    cpg_ONT, ONT_normal, ONT_tumor, cpg_bisulf, normal_bisulf, tumor_bisulf, out_fname = sys.argv[1:8]

    # Load data from files
    bisulfite=load_data(cpg_bisulf)
    bisulfite_n= load_data(normal_bisulf)
    #print(bisulfite_n)
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

    ONT_common_sites= ONT_bs_overlap.intersection(ONT_set)
    BS_common_sites= ONT_bs_overlap.intersection(bisulfite_set)


    #print(ONT_bs_overlap)

    #how many are only bisulfite reads
    print(f"bisulfite only reads: {len(bisulfite_set)-len(ONT_bs_overlap)}")
    #BS_only=len(bisulfite_set)-len(ONT_bs_overlap)
    #print(BS_only)

    #how many are ONT reads
    print(f"ONT only reads: {len(ONT_set)-len(ONT_bs_overlap)}")
    #ONT_only=len(ONT_set)-len(ONT_bs_overlap)



    #how many reads are shared between ONT and bisulfite reads
    print(f"percentage of shared reads: {len(ONT_bs_overlap)/len(ONT_bs_together)*100} %")


    fig, ax = plt.subplots(3, 1)
    #figsize=(5,15)

    #set the labels and the title =================================
    coverage_value_ont = numpy.array([x[1] for x in ONT.values()])
    ax[0].hist(coverage_value_ont, bins=numpy.amax(coverage_value_ont), color='blue', label = "ONT")
    coverage_value_bisulfite = numpy.array([x[1] for x in bisulfite.values()])
    ax[0].hist(coverage_value_bisulfite, bins=numpy.amax(coverage_value_bisulfite), color='orange', label = "bisulfite")
    ax[0].set_xlim(0,100)
    # label axes, change color, add title ==========================================
    
    #3C===========================================================================================================================================
    #position in overlap, find in ONT or BS the [0] value at the position
    ONT_overlap_sites=numpy.array([ONT[x][0] for x in ONT_bs_overlap])
    BS_overlap_sites=numpy.array([bisulfite[x][0] for x in ONT_bs_overlap])
    np_hist_overlap, x_edges, y_edges=numpy.histogram2d(ONT_overlap_sites, BS_overlap_sites, bins=100)
    ax[1].imshow(numpy.log10(np_hist_overlap+1))
    pearsons=numpy.corrcoef(ONT_overlap_sites, BS_overlap_sites)[0,1]
    title = f'Relationship of Methylation Calls Between ONT and BS Sequencing with Pearson\'s Correlation: {pearsons:.2f}'
    ax[1].set_title(title)
    #add in title with r value and axis titles 

    #plt.show()
    #0,1 in tuple
    
    #calc R pearson
    #google: include variable as string in a title to include the pearson's coeff in the title 




    #3D================================================================================================================

    n_t_change_BS_sites = {}
    # tumor_BS=[]
    # normal_BS=[]
    for key,value in bisulfite_n.items():
        #print(key,value)
        if key not in bisulfite_t:
            continue
        # else:
        #     tumor_BS.append(tumor_BS[value])
        #calculate change between normal tumor 
        change_BS = bisulfite_t[key][0]-value[0] #grabbing the 0th value associated with the key in BS tumor and subtracting 
        if change_BS != 0:
            n_t_change_BS_sites[key]=change_BS #assign at the key position the value of change in normal and tumor methylation


    n_t_change_ONT_sites = {}
    for key,value in ONT_n.items():
        #print(key,value)
        if key not in ONT_t:
            continue
        #calculate change between normal tumor 
        change_ONT = ONT_t[key][0]-value[0]
        if change_ONT != 0:
            n_t_change_ONT_sites[key]=change_ONT
    #print(n_t_change_BS_sites.values())



    common_change_sites=[]
    for key,value in n_t_change_BS_sites.items():
        if key in n_t_change_ONT_sites:
            common_change_sites.append((value, n_t_change_ONT_sites[key]))
    #common chnage sites now includes tuple with each common key where value is info from BS and ONT
    common_change_sites=numpy.array(common_change_sites)
    #print(common_change_sites)



    #plotting the changes - ADD THE TITLE AND AXIS LABELS ==============================
    ax[2].violinplot((list(n_t_change_ONT_sites.values()), list(n_t_change_BS_sites.values())))
    #ax[2].set_xticks([1,2])
    #ax[2].set_xticklabels("ONT", "Bisulfite")
    pearsons_2=numpy.corrcoef(common_change_sites[:,0], common_change_sites[:,1] )[0,1]
    title2 = f'Distribution of ONT vs BS Methylation Call Changes with Pearson\'s Correlation: {pearsons_2:.2f}'
    ax[2].set_title(title2)
 
   






    plt.tight_layout() 
    plt.show()









def load_data(fname):
    data={}
    for line in open(fname):
        line=line.rstrip().split("\t")   #\t split on tab
        key=(line[0], line[1])
        data[key]=(float(line[3]), int(line[4]))
    return(data)  







main()

































#TRASH 

 # list_ONT= list(ONT_set)
    # #print(list_ONT[:10])
    # flattened_list_ONT=[]
    # for i in list_ONT:
    #     for j in i:
    #         if j != "chr2":
    #             j=int(j)
    #             #print(type(j))
    #             flattened_list_ONT.append(j)
    # #print(flattened_list_ONT)
    # list_BS= list(bisulfite_set)
    # #print(list_BS[:10])
    # flattened_list_BS=[]
    # for i in list_BS:
    #     for j in i:
    #         if j != "chr2":
    #             j=int(j)
    #             #print(type(j))
    #             flattened_list_BS.append(j)
    # #print(flattened_list_BS[:10])




