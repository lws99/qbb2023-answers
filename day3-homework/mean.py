#!/usr/bin/env python

import sys

#EXERCISE 1: 1.) Without using any external libraries (such as numpy) write a function that takes a list 
#(of any length) of integers as input and returns the mean (i.e., average). 
#Write a script called mean.py where you create the list of integers, compute the mean using your function, and print it.



def mean_of_integers(list):
    sum_of_nums=0
    for number in list:
        sum_of_nums += number # adds the numbers together from list as they are added to sum_of_numbers
    average_of_list=sum_of_nums/(len(list)) #don't have this line and return in for loop or it will be repeated with every for loop run
    return average_of_list


list_of_integers=[2,2,2,2,2,2,2,2,2,2,2]
print(mean_of_integers(list_of_integers))






