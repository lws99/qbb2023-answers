#!/usr/bin/env python

import sys


#3.) Write a function that takes a patient row index (integer) and data file name (string) as 
#input and returns the mean inflammation level across the 40 days (float) for that given patient. 
#Embed this function in a script called mean_inflammation.py that defines a patient row index as a variable, executes the function, and prints the output.



def mean_inflammation(patient_index, data_file_name):
    dataset=open(data_file_name).readlines()  #this opens the user's
    patient_of_interest=dataset[patient_index]
    patient_of_interest=patient_of_interest.strip()
    patient_of_interest=patient_of_interest.split(",")
    day_list=[] 
    for day in patient_of_interest:
        day=int(day)
        day_list.append(day)
    return sum(day_list)/(len(day_list))

print(mean_inflammation(12, 'inflammation-01.csv'))     












