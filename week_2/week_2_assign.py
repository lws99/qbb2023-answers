#!/usr/bin/env python


    #Step 1.1

#In your README.md for this assignment, answer the following question (show your work):
#How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?

#Do 3 million divided by 100 *GOAL IS TO SIMULATE READ COVERAGE (runnning simulation at least 3 times --> write a function)


import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as stats
import statsmodels.formula.api as smf
import statsmodels.api as sm


#=====================================================================================================================================

    #STEP 1.2

def simulate_coverage(coverage, genome_len, read_len, figname):
    #simulate the genome reads as a numpy array with zeros (zeros until you simulate the reads)
    coverage_arr=np.zeros(genome_len)
    #print(coverage_arr)
    #figure out how many reads you need to generate as an integer --> no floats 
    num_reads=int(coverage*genome_len/read_len)
    print("number of reads required is:", (num_reads))
    
    #randomly pick a start position in the genome (using numpy function that randomly picks an integer between 2 numbers) and randomly increment through the rest of the genome 
    low=0 #(start position in genome --> position you can start from in genome is 0)
    high= genome_len - read_len #end position, the highest you can start from is the read length (you can't start from less than the read length at the end of the genome because then the reading would run off of the genome length)
    start_position=np.random.randint(low=0, high=high+1, size=num_reads) #low is inclusive and high is not inclusive (whhy you do +1 so that it does include the number high is), size tells numpy how many random integers you want
    print("number start position is:", len(start_position)) #each start position represents a read 

    #loop through the start positions so you have this coverage tracking for every single read --> increment start position coverage and all following reads coverage in the read len after the start postition by one
    for start in start_position:
        coverage_arr[start: start+read_len] +=1 #you want to increment not just the single start postion, but all of the following positions in the read length after the start
        

        #make histogram
    x=np.arange(0, max(coverage_arr)+1) #make a range of numbers from 0 to max seen in simluation (+1 because max is exclusive) --> max number of coverages seen in array
    #this moves the histogram over by 0.5 (just an asthetic choice, doesn't really matter)
    bins=x-0.5
    sim_0cov=genome_len-np.count_nonzero(coverage_arr) #how many zero coverage bases are present in the coverage array (what number of bases are not covered in genome)
    sim_0cov_percent=   100* (sim_0cov/genome_len)                   #what percentage of bases in genome are 0 (not covered)
    print(f'In this simulation, there are {sim_0cov} bases with a zero read')
    print(f'This is {sim_0cov_percent}% of the genome') #YOU CAN INCREASE PERCENTAGE OF GENOME COVERED BY INCREASING THE NUMBER OF READS 
    
    #get the poisson distribution 
    y_poisson=stats.poisson.pmf(x , mu=coverage)  *genome_len #pmf is probability a distribution will give you a particular number (mean of poisson distruction is lambda, but is mu in scipy.stats). Mu is the coverage that you are targeting in this case it's 3. We want to know what the probaility of getting a certain number of reads is as well (how many bases have the certain number of coverage --> this is x in this case because x is the max number of coverages seen in the array).  #find total area under the histogram (the area under the histogram is equal to the total genome length, in this case it is 1 million) aka genome len
    #the numbers in the positions of the prints tells you the probability of getting each number in the distribution (area under the poisson distribution is ~1, but not exactly 1 --> gives you probabilites)
    print(y_poisson)

    #get normal distribution
    y_normal=stats.norm.pdf(x, loc= coverage, scale = np.sqrt(coverage))*genome_len #Loc is the mean. you have to tell the normal distribution the mean and SD. mean is coverage and scale (SD) is the square root of the mean. multiply by genome_len
    print(y_normal) #these are proximations of probabilities (converts liklihoods to proximation of probablities )


    #make a histogram
    fig, ax = plt.subplots()
    ax.hist(coverage_arr, bins=bins, label="Simulation") #bins makes each coverage its own bin
    ax.plot(x, y_poisson, label="Poisson") #the poisson distribution tells us the expected coverage based on math of variables. Here the poisson matches pretty well with the distribution which means that the simulation is pretty good. 
    ax.plot(x, y_normal, label="normal") # this fits worse than the poisson distribution, but as coverage is increased the fit of the normal distribution (and Poisson) will get better --> bigger the mean the better the aproximation is 
    #x is number of reads and y is frequency of the bp
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Frequency (bp)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figname)
    


   # plt.show()


#NORMAL DISTRIBUTION: the PROBABILITY of any one event happening in a single distribution cannot be found (THIS CAN BE FOUND IF THE PROBABILITY OF FINDING SOMETHING IS IN A RANGE --> IS AREA UNDER CURVE IN RANGE) but LIKELIHOOD of an event happening can be found 


simulate_coverage(3, 1_000_000, 100, "ex1_3x_cov.png")

#how many reads are neccessary to do the above?
#30,000

#======================================================================================================================================================================
    #Step 1.4

#Now, repeat the analysis with 10x coverage:

#Simulate using the appropriate number of reads
#Make a histogram. Overlay a Poisson distribution with lambda = 10. Overlay a Normal distribution with mean = 10 and std. dev. = 3.16. 
#Upload this plot as ex1_10x_cov.png in your submission directory.
#In your README.md, answer the following questions:
#In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
#How well does this match Poisson expectations? How well does the normal distribution fit the data?


simulate_coverage(10, 1_000_000, 100, "ex1_10x_cov.png")


#======================================================================================================================================================================

    #Step 1.5

#Now, repeat the analysis with 30x coverage:

#Simulate using the appropriate number of reads
#Make a histogram. Overlay a Poisson distribution with lambda = 30. Overlay a Normal distribution with mean = 30 and std. dev. = 5.47. 
#Upload this plot as ex1_30x_cov.png in your submission directory.
#In your README.md, answer the following questions:
#In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
#How well does this match Poisson expectations? How well does the normal distribution fit the data?



simulate_coverage(30, 1_000_000, 100, "ex1_10x_cov.png")




#======================================================================================================================================================================
#======================================================================================================================================================================
                                                #EXERCISE 2  ####I AM STILL WORKING ON THIS AND KNOW THAT THIS IS NOT FULLY CORRECT!####

#======================================================================================================================================================================
#======================================================================================================================================================================



    #Step 2.1
#Next, you’re going to generate your own de Bruijn graph using a provided set of reads. Copy the list of reads below into your code:
#Write code to find all of the edges in the de Bruijn graph corresponding to the provided reads using k = 3 (assume all reads are from the forward strand, 
#no sequencing errors, complete coverage of the genome). Each edge should be of the format ATT -> TTC. Write all edges to a file, with each edge as its own line 
#in the file.


reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']


graph=set()


k=3
for read in reads:
  for i in range(len(read) - k):
     kmer1 = read[i: i+k]
     #print(kmer1)
     kmer2 = read[i+1: i+1+k]
     #print(kmer2)
     graph.add(kmer1)
     graph.add(kmer2)
     

#print(graph)


for edge in graph:
   print(edge)















#commonly oversample 30X to be confident that you are covering enough of the genome


















