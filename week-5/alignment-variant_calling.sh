# Step 1.1: Index the sacCer3 genome

# You’ll be using a tool called bwa (BWA manual) to perform alignment. Before you can align your sequencing reads, 
# bwa needs you to index the sacCer3 genome. Without getting into the the nitty gritty, this essentially means creating a table of contents (and index) 
# that bwa can use to quickly find matches between your reads and the reference genome.

# Using bwa index, create an index for the sacCer3.fa reference.

#==================DONE===================================================




# Step 1.2: Align your reads to the reference

# Now that you’ve indexed the reference, you can align your reads to the reference using bwa mem. Because you’ll want to run this 
# step the same way on all 10 strains, it makes sense to do this step in a bash for loop. Consider this resource for how to write a 
# for loop in a bash script: bash for loops walkthrough.

# Create a bash for loop that loops through each of the 10 samples. For each sample, use bwa mem to align the reads to the reference.

# IT IS VERY IMPORTANT that you assign each sample a read group during this process, so that individual samples can be distinguished later 
# in Step 2.1. You can do this with the (somewhat cryptic) -R flag, which you use to add a line to the header of each output alignment file. 
# An example of a header line you can add with the -R flag is "@RG\tID:Sample1\tSM:Sample1". You can replace “Sample1” here with the appropriate 
# sample name for each of your yeast strains.

# Perhaps consider the -t and -o flags as well.




for sample in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
do
    echo "Aligning sample:" ${sample}
    bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}" \
      sacCer3.fa \
      ${sample}.fastq > ${sample}.sam
      samtools sort ${sample}.sam -o ${sample}.bam
      samtools index ${sample}.bam
# done

#=====================================DONE======================================================================================================


# Step 1.3: Format and index your alignments

# Now that you’ve aligned your reads to the reference, you should have 10 .sam files, one for each sample. 
# These files contain all of the alignments for each yeast strain. You can see how they’re organized with less -S, and you can read more about the SAM format here.

# These files contain all of the information you need for variant calling, but before you can do that, they’ll need to be sorted and indexed 
# (similar to how you indexed the reference in Step 1.1. For both of these tasks you can use the samtools program (manual here, or you can just run samtools help).

# First, sort each of your .sam files using samtools sort. You can do this in a new for loop in your bash script or, even better, in the same for 
# loop you used for alignment. You’ll want to output these sorted files as .bam files, which contain the same information as the .sam file but are compressed.

# Perhaps consider the -O and -o flags when running samtools sort. --> allows you to turn the .sam into .bam

# Next, create an index for each of the resulting sorted .bam files using samtools index. As before, you can do this in a new for loop in your bash 
# script or in the same for loop as the previous two steps.

# At the end of this step, you should have 10 sorted .bam files and their corresponding .bam.bai indices.

#=====================================DONE=======================================================================================================




# Step 2.1: Call variants

# For variant calling, you’ll be using a tool called freebayes (manual here, or you can just run freebayes --help).

# Use freebayes to identify genetic variants in all of your yeast strains concurrently (i.e. you should only be running freebayes once will all samples, 
#     not for each sample separately). It will output results in Variant Call Format (.vcf). You can read more about the VCF format here.

# You should consider using the -f, --genotype-qualities, and -p flags. You might like the -L flag as well.

# NOTE: For TAs, running this step took nearly 15 minutes. We expect this step will take a similar amount of time for you, and your 
# computer might make a lot of noise.



freebayes -f sacCer3.fa A01_09.bam A01_11.bam A01_23.bam A01_24.bam A01_27.bam A01_31.bam A01_35.bam A01_39.bam A01_62.bam A01_63.bam --genotype-qualities > A01.vcf




#=======================================DONE=====================================================================================================


# Step 2.2: Filter variants based on site quality

# Sometimes, freebayes will call variants that may not be “real”. Luckily, variant callers generally report, for each variant, a “site quality”, that 
# describes the probability of that site truly being polymorphic. You can use that quality score to filter out low quality (low confidence) variants.

# Note that the “site quality” is not the same as the “genotype quality”, which describes the probability that a single sample’s genotype at some variant 
# is correct. Both qualities, however, are generally reported on the “Phred” scale (more here).

# You can filter out low quality variants using the vcffilter tool (documentation here).

# Filter your VCF using vcffilter so that you only keep variants whose estimated probability of being polymorphic is greater than 0.99. 
# You should consider how to do this with the -f flag. Output your filtered variant calls to a new VCF file.


#filter qual>20 to find the best quality site reads 



vcffilter -f "QUAL > 20" A01.vcf > filtered_A01.vcf




#==================================DONE===========================================================================================================


# Step 2.3: Decompose complex haplotypes

# Sometimes, especially in regions with complex alignments between samples, you can run into cases where there are multiple possible alternative 
# alleles at the same position. There’s nothing inherently “wrong” with this, but it does often complicate downstream analyses.

# Luckily, you can use the vcfallelicprimitives tool (documentation here) to decompose these more complex haplotypes into more manageable biallelic variants.

# Use vcfallelicprimitives to decompose complex haplotypes in your filtered VCF from Step 2.2. We suggest using the -k and -g flags to keep annotations 
# for the variant sites and sample genotypes in your VCF. Output to a new VCF file.



vcfallelicprimitives -k -g filtered_A01.vcf > decomposed_filtered_A01.vcf

less -S decomposed_filtered_A01.vcf 

#=====================================DONE============================================================================================




# Step 2.4: Annotate variants

# Now that you’ve got these high-quality and nicely behaving variant calls, you want to know what impact these variants might have. Obviously, this is a huge question, 
# and is the basis of all of genetics, but we can get a basic idea of their functional impact on nearby genes (e.g. are they missense, nonsense, etc.) using the snpEff 
# tool (the online documentation isn’t great, but running snpEff ann --help should provide some useful information).

# If it wasn’t obvious, snpEff requires prior annotations (e.g. gene annotations) to work. Have snpeff download its database of Saccharomyces cerevisiae annotations 
# using the following command (we told you the NCBI ID would be relevant):

# snpEff download R64-1-1.105
# Finally, use snpEff ann to annotate your VCF with the predicted functional effects that these genetic variants may have. Output to a new (and final) VCF.

# For submission purposes, use head to grab just the first 100 lines of your final VCF and store this in a new VCF. You will submit this “sample” VCF along with the 
# rest of your assignment. YOU SHOULD NOT SUBMIT ANY OTHER VCFS; THEY ARE TOO BIG. Depending on how your .gitignore is set up, you may need to do 
#     git add --force <yoursamplevcf.vcf>.




#I WILL NEED TO DO THIS WHEN I PUSH TO GITHUB!!!!               git add --force <yoursamplevcf.vcf>



#Genome variants in yeast:
snpEff download R64-1-1.105


#snpEff ann help
snpEff ann R64-1-1.105 decomposed_filtered_A01.vcf > annotated_decomp_filt_A01.vcf


#command line              head -100 decomposed_filtered_A01.vcf |

#head -100 annotated_decomp_filt_A01.vcf > 100_A01.vcf




