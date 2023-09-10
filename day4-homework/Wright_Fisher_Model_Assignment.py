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


#Creating a function that will return allele frequencies for each generation of a population of at least 100
def W_F_model(starting_allele_freq, pop_size=100):
    #find the allele freq of the next pop
    allele_freq=[]
    while starting_allele_freq < 1 and starting_allele_freq > 0:
        success=np.random.binomial(2*pop_size, starting_allele_freq) #n is the number of chromosomes
        starting_allele_freq= success/(2*pop_size)
        allele_freq.append(starting_allele_freq)
    return allele_freq


#finding the allele frequencies for an allele with  starting frequency of 0.8 in a population of 100 using W_F_model
allele_frequencies_test=W_F_model(0.5)
#print(allele_frequencies_test)

#finding the number of generations it will take for the allele to reach fixation
print(len(allele_frequencies_test))




# Plot data


#preparing the variables for plotting
generation=range(len(allele_frequencies_test))
frequency=allele_frequencies_test



#plotting time to allele fixation for one model iteration
fig0, ax = plt.subplots()
ax.set_title( "Time to allele fixation" )
ax.set_xlabel("generation")
ax.set_ylabel("frequency")
ax.plot(generation, frequency, c = "red")







#EXERCISE 2

#plotting allele frequency trajectories across generations for 30 model iterations
fig1, ax=plt.subplots()
for i in range(30):
    allele_frequencies_test=W_F_model(0.5)
    generation=range(len(allele_frequencies_test))
    frequency=allele_frequencies_test
    ax.set_title( "Time to allele fixation" )
    ax.set_xlabel("generation")
    ax.set_ylabel("frequency")
    ax.plot(generation, frequency)




#making a histogram of times to fixation
fig2, ax=plt.subplots()
for i in range(1000):
    allele_frequencies_test=W_F_model(0.5)
    frequency=allele_frequencies_test
    ax.set_title( "Distribution of time to allele fixation" )
    ax.set_ylabel("frequency")
ax.hist(frequency)












#EXERCISE 3

#Pick at least five population sizes greater than or equal to 50. For each population size, run the model at least 50 times and find the average time to fixation. 
#Keep your allele frequency constant for all runs. Create a scatter or line plot of population size vs. average time to fixation.

#affect of population size on allele fixation


#create a list of different population sizes
size_of_population=[50,60,70,80,90]
#create an empty list where the mean time to fixation for each population size will be stored 
mean_final=[]
#running the model 50 times on each population size
for p in size_of_population:
    for i in range(50):
        generations_3=(len(W_F_model(0.5, p)))
#taking the mean of time to fixation for each population
    mean=np.mean(generations_3)
    mean_final.append(mean)


#list of the mean time to fixation for each population in size_of_population
print(mean_final)



#plotting the data size vs average time to fixation
fig3, ax=plt.subplots()

population=size_of_population
time_fixation=mean_final
ax.set_title( "population size vs mean time to fixation")
ax.set_xlabel("population size")
ax.set_ylabel("mean time to fixation")
ax.scatter(population, time_fixation)



#This time, pick a population size and vary the allele frequency. Run at least 10 trials for each allele frequency.
#affect of allele frequency on fixation


#create a list of different frequencies
diff_starting_freqs=[0.2,0.4,0.6,0.8]
#create an empty list where the mean time to fixation for each starting frequency will be stored 
mean_final_freqs_fixation=[]
#running the model 10 times on each starting frequency
for starting_freq in diff_starting_freqs:
    for i in range(10):
        frequencies_3=(len(W_F_model(starting_freq, 500)))
#taking the mean of time to fixation for each starting allele frequency
    mean_time_to_fixation_for_frequencies=np.mean(frequencies_3)
    mean_final_freqs_fixation.append(mean_time_to_fixation_for_frequencies)
print(mean_final_freqs_fixation)



#plotting the starting frequencies vs mean time to fixation
fig4, ax=plt.subplots()

changing_starting_freqs=diff_starting_freqs
time_to_fixation=mean_final_freqs_fixation
ax.set_title( "starting allele frequency vs mean time to fixation")
ax.set_xlabel("starting allele frequency")
ax.set_ylabel("mean time to fixation")
ax.scatter(changing_starting_freqs, time_to_fixation)


plt.show()


plt.close( fig0 )
plt.close( fig1 )
plt.close( fig2 )
plt.close(fig3)
plt.close(fig4)


