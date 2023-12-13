Link to any outside resources in this README



Q1: Are the majority of the CpG dinucleotides methylated or unmethylated?
The majority are methylated.




Q2: How does using nanopore for methylation calling differ from bisulfite sequencing in terms of coverage? Which method appears better and why?
Nanopore has less coverage than bisulfite sequencing. In terms of coverage, bisulfite sequencing is best. 

3A #2:
Number of sites in bismark: 132466
Number of sites in nanopore: 53346
Percentage of shared sites:95.85468138657977 %



Question 2:  
Nanopore sequencing has deeper coverage than bisulfite sequencing. Overall, I think Nanopore sequencing is the best because it covers the genome with deeper coverage than bisulfite sequencing, giving you more information about methylation across the entire genome. 

Question 3:
The ability to detect methylation changes between the two approaches is essentially the exact same, especially because Pearson's Coefficient is 1 (perfect correlation). Tumorigenesis tends to increase methylation at CpG sites. 


Question 5:
In the tumor landscape, there is overall more DNMT3A methylation than in the normal landscape. This is shown in the bisulfite sequencing having more peaks in the tumor sample, and the ONT sample appearing to have broadly more methylation in the tumor sample than the normal sample. This increase in DNMT3A methylation likely overall decreases de novo methylation in tumor samples. 


Question 6:

Gene imprinting means that a gene is differentially expressed based on the parent that a gene is inherited from. When I try to phase the reads the program tells me there are no variants in the selected range. 

Question 7:
Variants are required for phasing, which means that more information (longer read sequencing) is required.


Quation 8:
No, not any set of reads can be phased because not every set contains variants. 
