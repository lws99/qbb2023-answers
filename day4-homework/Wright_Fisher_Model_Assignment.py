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








#EXERCISE 1

#PART1: Creating a function that will return allele frequencies for each generation of a population of at least 100
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



#PART2: finding the number of generations it will take for the allele to reach fixation
print(len(allele_frequencies_test))



#PART3: Plotting allele trajectory --> plt.show() is at the end of my code, so when you run the entire thing it will show all of the graphs

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

#PART1: plotting allele frequency trajectories across generations for 30 model iterations --> plt.show() is at the end of my code, so when you run the entire thing it will show all of the graphs
fig1, ax=plt.subplots()
for i in range(30):
    allele_frequencies_test=W_F_model(0.5)
    generation=range(len(allele_frequencies_test))
    frequency=allele_frequencies_test
    ax.set_title( "Time to allele fixation" )
    ax.set_xlabel("generation")
    ax.set_ylabel("frequency")
    ax.plot(generation, frequency)




#PART2: making a histogram of times to fixation --> plt.show() is at the end of my code, so when you run the entire thing it will show all of the graphs
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


#PART1: plotting the population size vs average time to fixation -->plt.show() is at the end of my code, so when you run the entire thing it will show all of the graphs
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



#PART2: plotting the starting frequencies vs mean time to fixation -->plt.show() is at the end of my code, so when you run the entire thing it will show all of the graphs
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





#BASIC EXERCISE: INTERPRETING DATA
#For any plot, explain the results you see. What might be contributing to it? What does it mean biologically?
#For any assumption in the Wright-Fisher model, how might changing that assumption affect the result? How might nature and biology violate these assumptions?


#ANSWER 1 EXPLAINING THE RESULTS OF FIG3 (exercise 3 part 1)
#For my plot of population size vs time to fixation (fig3), what I expected to see was that larger populations should take more time to fix an allele than 
#smaller populations. My prediction is not always true when I run the Wright-Fisher Model in exercise 3, 
#likely due to the fact that the population sizes I chose are not that different from each other, so their resepective times to fixation may not be that different. 
#If I chose a much higher population size, the behavior of the model would likely be more predictable, and time to fixation for an allele would 
#essentially always take longer than for the smaller populations, like the ones I used. This is because larger populations experience less overall allele variation than 
#smaller populations. 




#ANSWER 2 ON THE ASSUMPTIONS OF THE MODEL- FIXED POPULATION
#If the assumption that there is constant population size between generations was not true within the model, then that would mean that the fluctuation of alleles 
#between each generation would not actually be random. Rather, allele fluctuation would likely be based on population size in each generation, which was demonstrated in 
# exercise 3 part 1 where changing population sizes had an effect on time to fixation. Nature and biology violates this assumption because every generation
#does not have a constant population size. 



#ANSWER 3 ON THE ASSUMPTIONS OF THE MODEL- NO SELECTION
#If the assumption that there is no natural selection was not true within the model, then alleles would not be able to be seen as evolutionarily neutral.
#This is because some allele variations are deadly or significantly worse for survival than other allele variations that are passed through generations.
#Thus, you would not be able to broadly understand an allele's fixation because of the high levels of outside pressure from selection agents, especially because
#selection is very dynamic and individualized depending on the population you are looking at. Natural selection also affects other variables that are fixed
# in the Wright-Fisher Model, like population, which further convolutes the fixation of an allele with selection present. 
#In nature/biology, not every allele is created equal in terms of the outcome for the organism, but the Wright-Fisher Model assumes that this is true. 


