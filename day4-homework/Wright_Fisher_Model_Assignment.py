#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#BERNOULLI TRIAL  : decribing random events
    #two outcoomes: success or failure, fixed probability over time (p)
    #more trials = tighter distributions

#numpy.random.binomial(n,p) performs Bernoulli trial

#pseudo code:

#get a starting frequency and population size
#input parameters for function
#use a while loop because you don't know how long the allele will take to reach fixation

#make a list to store our allele frequencies

#While ur allele frequency is between 0 and 1:
#   Get the new allele frequency for the next generation
#   by drawing from the binomial distribution  --> this will give you number of successes
#   (covert number of successes into a frequency)

#   Store our allele frequency in the allele frequency list



#Return a list of allele frequency at each time point
#The number of generations to fixation 
#is the length of your list

#Your function should run until one of the two alleles reaches fixation (i.e., your allele frequency hits 0 or 1). #looking at two different alleles in a population and trying to figure out which one goes to fixation
#np.random.binomial(n,p)



def W_F_model(starting_allele_freq, pop_size=100):
    #find the allele freq of the next pop
    allele_freq=[]
    while starting_allele_freq < 1 and starting_allele_freq > 0:
        success=np.random.binomial(2*pop_size, starting_allele_freq) #n is the number of chromosomes
        starting_allele_freq= success/(2*pop_size)
        allele_freq.append(starting_allele_freq)
    return allele_freq


#finidng the allele frequencies for an allele with  starting frequency of 0.8 in a population of 100 using W_F_model
allele_frequencies_test=W_F_model(0.5)
#print(allele_frequencies_test)

#finding the number of generations it will take for the allele to reach fixation
print(len(allele_frequencies_test))




# Plot data
#preparing the variables for plotting
generation=range(len(allele_frequencies_test))

frequency=allele_frequencies_test



fig, ax = plt.subplots()
ax.set_title( "Time to allele fixation" )
ax.set_xlabel("generation")
ax.set_ylabel("frequency")
ax.plot(generation, frequency, c = "red")


plt.show()
plt.close( fig )




#EXERCISE 2




