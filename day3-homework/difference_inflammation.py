#!/usr/bin/env python


#Write a function that takes two patient row indices (integers) and data file name (string) as input and returns a list of the difference between their 
#inflammation levels on each of the 40 days (floats). 
#Embed this function in a script called difference_inflammation.py that defines patient row indices as variables, executes the function, and prints the output.


def patient_diffs(patient_row_1, patient_row_2, data_file_name):
    dataset=open(data_file_name).readlines()

    patient_1=dataset[patient_row_1]
    patient_1=patient_1.rstrip()
    patient_1=patient_1.split(",")

    patient_2=dataset[patient_row_2]
    patient_2=patient_2.rstrip()
    patient_2=patient_2.split(",")

    patient_1_int=[]
    for day in patient_1:
        day=int(day)
        patient_1_int.append(day)
    patient_2_int=[]
    for day in patient_2:
        day=int(day)
        patient_2_int.append(day)

    difference_between_inflammation=[]


    index=0
    for value in range(len(patient_1_int)):
        difference=(patient_1_int[index])-(patient_2_int[index])
        difference_between_inflammation.append(difference)
        index+=1
        


    return difference_between_inflammation




print(patient_diffs(0,1,'inflammation-01.csv'))
print(patient_diffs(5,40,'inflammation-01.csv'))


