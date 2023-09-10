#!/usr/bin/env python

import numpy

dataset=open("inflammation-01.csv", "r")
#print(dataset)
patients=dataset.readlines()
#print(patients)




##Exercise 1
patients_alone=[]
for patient in patients:
	patient=patient.rstrip()
	#print(patient)
	patient_list=patient.split(',')
	integer_patients=[] #starts as a list (partient_list) inside of a list (integer_patients) and we need to change it to in integer (of patient_list) inside of list integer_patients everytime
	for i in patient_list:
		i = int(i)
		integer_patients.append(i)
	patients_alone.append(integer_patients)
	
#indexing the total patient list (as int) for the fifth patient
fifth_patient = patients_alone[4]

#indexing the fifth patient's list of int for the flare-ups on day 1, 10, and the last day
first_day_flare_up=fifth_patient[0]
tenth_day_flare_up=fifth_patient[9]
last_day_flare_up=fifth_patient[-1]


#FINAL ANSWER: Fifth patient's flare ups on day 1, 10, and the last day
print(first_day_flare_up)
print(tenth_day_flare_up)
print(last_day_flare_up)









##Exercise 2 
##For each patient, calculate the average number of flare-ups per day. Print the average values for the first 10 patients.
##These are the row averages - for example, patient 1 has 5.45 flare-ups per day on average; patient 2 has 5.425 flare-ups per day on average.

individual_flares=[]
for patient in patients_alone:
	individual_flares.append(numpy.mean(patient))

#average individual flare-ups per day for each of the individual patients
print(individual_flares)

#average individual flare-ups per day for patients 1-10
print(individual_flares[1:10])
	








##Exercise 3 Using the average flare-ups per day calculated in part 2, print the highest and lowest average number of flare-ups per day.
#FINAL ANSWER: min and max of the average number of flare-ups per day
print(numpy.max(individual_flares))
print(numpy.min(individual_flares))









##Exercise 4 For each day, print the difference in number of flare-ups between patients 1 and 5.
first_patient = patients_alone[0]
fifth_patient


differences_between_one_and_five=[]
index=0
for day in range(len(first_patient)):
	difference= (first_patient[index])-(fifth_patient[index])
	differences_between_one_and_five.append(difference)
	index +=1

#FINAL ANSWER:  printing the difference in number of flare-ups, per day, between patients 1 and 5 
print(differences_between_one_and_five)
