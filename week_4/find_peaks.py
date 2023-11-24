#!/usr/bin/env python

import sys

from ChIP_seq_live_code import load_bedgraph, bin_array
import numpy
import scipy.stats
import matplotlib.pyplot as plt


#control sample and test sample
#want to normalize them to the same height (transform the the control peak to have the same amplitude as the test sample --> do it based on the coverage of each sample)
# number of reads that terminate at a positon on the chromosome --> use this to get coverage (add the number of reads to get total coverage of genome, divide the number of reads by the coverage)

#ChIP-seq data were mapped against the dm6 genome and were preprocessed to pull out the end location of each read, restricted to chromosome chr2R, partitioning them by strand, and creating pile-ups in bedgraph format. There is also a bed file containing all peak calls from chr2R, as your analysis will be restricted to a smaller region for processing time considerations.

#There are also two python scripts to serve as a framework for how to approach the assignment. They outline the general steps in comments and include a couple of functions you will need to accomplish each task.


#=======================================================================================================================================================================
                                            #Step 1: Calling peaks
#Using the estimated fragment size determined during the live coding portion of class and the bedgraph files for each sample and the control, 
#determine regions with p-values <= 1e-10 given the local background under the poisson distribution. 
#You will be using the find_peaks.py script as your framework for this portion. Your script should process one sample at a time and you will run it for each sample. 
#Restrict your analysis to the region chr2R:10,000,000-12,000,000. Your script should produce the following outputs:

#A wiggle file with base-pair resolution -log10-transformed p-values
#A bed file with all regions with p-values <= 1e-10

from ChIP_seq_live_code import load_bedgraph

def main():
    # Load file names and fragment width
    forward_fname_sample, reverse_fname_sample, forward_fname_control, reverse_fname_control, output_prefix, frag_width = sys.argv[1:7]
    frag_width=int(frag_width)


    #./find_peaks.py sample1.fwd.bg sample1.rev.bg control.fwd.bg control.rev.bg 198

    #frag_width is 198 which is the offset that we found in the live code (it is the size of the ChIP seq peak at your seq of interest)
    
    # Define what genomic region we want to analyze
    chrom = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart




    # Load the sample bedgraph data, reusing the function we already wrote
    sample_forward=load_bedgraph(forward_fname_sample, chrom, chromstart, chromend)
    #print(numpy.mean(sample_forward))
    sample_reverse=load_bedgraph(reverse_fname_sample, chrom, chromstart, chromend)



    #Combine tag densities, shifting by our previously found fragment width
    sample=numpy.zeros(chromlen,float)
    sample[:-frag_width//2]+=sample_reverse[frag_width//2:]
    sample[frag_width//2:]+=sample_forward[:-frag_width//2]
    #print(numpy.mean(sample))



    # Load the control bedgraph data, reusing the function we already wrote
    control_forward=load_bedgraph(forward_fname_control, chrom, chromstart, chromend)
    control_reverse=load_bedgraph(reverse_fname_control, chrom, chromstart, chromend)



    # Combine tag densities, shifting by our previously found fragment width
    control=numpy.zeros(chromlen, float)
    control[:]=control_reverse[:]
    control+=control_forward[:]
    #print(numpy.mean(control))
  
    

    # Adjust the control to have the same coverage as our sample --> ie scaling
    control*=numpy.sum(sample)/numpy.sum(control)  #normalizing the control data to the sample data
    #print(numpy.mean(control))

    # Create a background mean using our previous binning function and a 1K window (binsize)
    # Make sure to adjust to be the mean expected per base
    background_mean=bin_array(control, 1000)/1000 #1000 walks through each base one at a time to find mean in the 1000

    #background you compare to for every new sample



    # Find the mean tags/bp and make each background position the higher of the
    # the binned score and global background score

    #background_mean=numpy.max(background_mean, control)
     #value in background or control higher?
    

    for i in range(len(background_mean)):
        if background_mean[i] < control[i]:
            background_mean[i] = control[i]


        



    # Score the sample using a binsize that is twice our fragment size
    # We can reuse the binning function we already wrote

    sample_score=bin_array(sample, 2*(frag_width))




    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background

    p_val= scipy.stats.poisson.cdf(sample_score, 2*background_mean*frag_width)
    print((("position pval", numpy.mean(p_val))))





    # Transform the p-values into -log10
    # You will also need to set a minimum pvalue so you doen't get a divide by
    # zero error. I suggest using 1e-250

    for i in range(len(p_val)):
        if p_val[i]>1e-250:
            p_val[i] = -numpy.log10(p_val[i])

    print(("log10 pval", numpy.mean(p_val)))




    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"

    write_wiggle(p_val, chrom, chromstart, f"{output_prefix}.wig")



    write_bed(p_val, chrom, chromstart, chromend, frag_width, f"{output_prefix}.bed") #take  output_prefix and use as a string, even though it is a variable I defined

def write_wiggle(pvalues, chrom, chromstart, fname):
    output = open(fname, 'w')
    print(f"fixedStep chrom={chrom} start={chromstart + 1} step=1 span=1",
          file=output)
    for i in pvalues:
        print(i, file=output)
    output.close()



def write_bed(scores, chrom, chromstart, chromend, width, fname):
    chromlen = chromend - chromstart
    output = open(fname, 'w')
    while numpy.amax(scores) >= 10:
        pos = numpy.argmax(scores)
        start = pos
        while start > 0 and scores[start - 1] >= 10:
            start -= 1
        end = pos
        while end < chromlen - 1 and scores[end + 1] >= 10:
            end += 1
        end = min(chromlen, end + width - 1)
        print(f"{chrom}\t{start + chromstart}\t{end + chromstart}", file=output)
        scores[start:end] = 0
    output.close()




if __name__ == "__main__":
    main()













# #Step 2: Intersecting peaks
# Because you have replicates, you can boost the confidence of your peak calls to those appearing in both samples. 
# To do this, you will be using bedtools, a program that performs genomic arithmetic. To get a list of peaks appearing in both samples, 
# use the following command:

# bedtools intersect -a sample1_peaks.bed -b sample2_peaks.bed > combined_peaks.bed
# Obviously you will need to change the file names if you use a different naming scheme for your peak calling output.

# Question: What fraction of peaks were retained when intersected compared to sample1? Sample2?

# sample1: 36%
# sample2: 36%



