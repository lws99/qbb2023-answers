#!/usr/bin/env python
import sys

#HOMEWORK 3 EXERCISE 2 
#2.) Use the function from part 1 to write a Python script called mean_from_file.py that computes and prints the mean of a series of integers from a data file 
#(e.g., “my_integers.txt”) where the data contain a set of integers, one per line. 
#[Optional: Write your code so that the user can specify the name of that file from the command line (hint: use sys.argv).]

#loading in the data from my file, using the command line
fname=sys.argv[1]
data=open(fname).readlines()


#cleaning the data from my dataset and converting the data type to integers
integer_data=[] #this needs to be outside of the for loop so that it does not get redefined every time that the for loop runs. If it was inside the for loop, then integer_data would just equal the last number to be read into the loop.
for line in data:
    stripped_line=line.rstrip()
    #split_line=stripped_line.rsplit() #the data does not need splitting here because we want all the lines to be read into one list
    data_stripped=stripped_line
    #for i in data: #this was unneeded because data_stripped is what we want to turn into integers
    i=int(data_stripped)
    integer_data.append(i)
    #print(integer_data)



#redefining the function from mean.py so that i can use it on integer_data from above
def mean_of_integers(list):
    sum_of_nums=0
    for number in list:
        sum_of_nums += number # adds the numbers together from list as they are added to sum_of_numbers
    average_of_list=sum_of_nums/(len(list)) #don't have this line and return in for loop or it will be repeated with every for loop run
    return average_of_list


#FINAL ANSWER: printing out the mean of my dataset that was loaded using the command line
print(mean_of_integers(integer_data))