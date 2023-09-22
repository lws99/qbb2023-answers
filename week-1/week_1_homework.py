#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.stats as sps
import statsmodels.formula.api as smf
import statsmodels.api as sm


#did a .gitignore with *.csv and the file names in the same directory as this .py -->(week-1) so that when i push this file to GitHub it will not push ANY of the .csv files



    #STEP 1.1
#You’ll start by exploring the data in aau1043_dnm.csv. First, load this data into a pandas dataframe.

df_offspring = pd.read_csv("aau1043_dnm.csv")
#print(df_offspring)
#column names are: Chr,Pos,Ref,Alt,Proband_id,Phase_combined,Crossover,Sanger



    #STEP 1.2
#You first want to count the number of paternally and maternally inherited DNMs (de novo mutations) in each proband. 
#Using this dataframe, create a dictionary where the keys are the proband IDs and the value associated with each key is a list of length 2, 
#where the first element in the list is the number of maternally inherited DNMs and the second element in the list is the number of paternally 
#inherited DNMs for that proband. You can ignore DNMs without a specified parent of origin.
#THE DNMs ARE HOW MANY TIMES EACH PROBAND ID SHOWS UP FOR EACH PARENT (proband IDs can show up multiple times)



#indexing the maternal and paternal data (boolean) and storing each in a list
maternal_proband_IDs=df_offspring.loc[:, "Phase_combined"] == "mother"
paternal_proband_IDs=df_offspring.loc[:,"Phase_combined"] =="father"

#print(maternal_proband_IDs)
#print(paternal_DNMs_rows)

#indexing all of the rows associated with each parent from Phase_combined
df_maternal=df_offspring.loc[maternal_proband_IDs, :]
#print(df_maternal)
df_paternal=df_offspring.loc[paternal_proband_IDs, :]
#print(df_paternal)



#make an empty dictionary with the key being Proband_id from df_offspring
offspring_dictionary={}
for proband_ID in df_offspring.loc[:, "Proband_id"]:
#setting the keys in the dictionary to proband ID and setting the values of the keys to 0,0
    offspring_dictionary[proband_ID] = [0, 0]


#add maternal Proband_id counts to the offspring_dictionary as the first value in the list 
for maternal_proband in df_maternal.loc[:, "Proband_id"]:
    offspring_dictionary[maternal_proband][0] +=1
    

#add father Proband_id counts to the offspring_dictionary as the second value in the list
for paternal_proband in df_paternal.loc[:, "Proband_id"]:
    offspring_dictionary[paternal_proband][1] +=1








    #STEP 1.3
#Use the following code snippet to convert this dictionary into a new pandas dataframe (this assumes your dictionary from step 1.2 is called deNovoCount):
#deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])
#Feel free to ask questions about how this code is working or, if you’re interested, you can try to figure it out yourself.

offspring_dictionary_DF=pd.DataFrame.from_dict(offspring_dictionary, orient="index", columns=["maternal_dnm", "paternal_dnm"])

#print(offspring_dictionary_DF)






    #STEP 1.4
#Now, load the data from aau1043_parental_age.csv into a new pandas dataframe.
#HINT: You will probably want to use the index_col argument with pd.read_csv(). It will make your life easier in the next step.


df_parental = pd.read_csv("aau1043_parental_age.csv", 
    index_col="Proband_id")
#print(df_parental)
#column names are: Proband_id,Father_age,Mother_age






    #STEP 1.5

#You now have two dataframes with complementary information. It would be nice to have all of this in one data structure. 
#Use the pd.concat() function (more here) to combine your dataframe from step 3 with the dataframe you just created in step 4 to create a new merged dataframe.
#HINT: You will want to specify the axis and join arguments with pd.concat()


concat_DF=pd.concat([df_parental, offspring_dictionary_DF], axis=1, join="inner")
#print(concat_DF)





                    #EXERCISE 2
#Using the merged dataframe from the previous section, you will be exploring the relationships between different features of the data. 
#The statsmodels package (more here) is an incredibly useful package for conducting statistical tests and running regressions. 
#As such, it is especially appropriate for the types of questions we’re interested in here. For this assignment, we’ll be using the formula api from 
#statsmodels (more here) to run some regressions between variables in our dataset. You can load this tool into Python with import statsmodels.formula.api as smf.





    #STEP 2.1
#First, you’re interested in exploring if there’s a relationship between the number of DNMs and parental age. Use matplotlib to plot the following. 
#All plots should be clearly labelled and easily interpretable.

#the count of maternal de novo mutations vs. maternal age (upload as ex2_a.png in your submission directory)
maternal_age=concat_DF.loc[:, "Mother_age"]
#print(maternal_age)
maternal_DNMs_alone=concat_DF.loc[:, "maternal_dnm"]
#print(maternal_DNMs)

fig1, ax = plt.subplots()
ax.set_title("Maternal age vs Maternal DNMs")
ax.set_xlabel("Maternal Age")
ax.set_ylabel("Maternal DNMs")


ax.scatter(maternal_age, maternal_DNMs_alone, c = "red")
#plt.show()


#the count of paternal de novo mutations vs. paternal age (upload as ex2_b.png in your submission directory)
paternal_age=concat_DF.loc[:, "Father_age"]
#print(paternal_age)
paternal_DNMs_alone=concat_DF.loc[:, "paternal_dnm"]
#print(paternal_DNMs)

fig2, ax = plt.subplots()
ax.set_title("Paternal age vs Paternal DNMs")
ax.set_xlabel("Paternal Age")
ax.set_ylabel("Paternal DNMs")


ax.scatter(paternal_age, paternal_DNMs_alone, c = "blue")
#plt.show()







    #STEP 2.2
#Now that you’ve visualized these relationships, you’re curious whether they’re statistically significant. 
#Perform ordinary least squares using the smf.ols() function to test for an association between maternal age and maternally inherited de novo mutations. 
#In your README.md for this assignment, answer the following questions:
#What is the “size” of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 2.1?
#Is this relationship significant? How do you know?

maternal_model=smf.ols(formula="maternal_dnm~1+Mother_age", data=concat_DF)
#running the ordiinary least squares and storing as a variable
maternal_results=maternal_model.fit()
#print(maternal_results.summary())
#print("maternal p-value is:", maternal_results.pvalues[1])







    #STEP 2.3
#As before, perform ordinary least squares using the smf.ols() function, but this time to test for an association between paternal age and paternally inherited 
#de novo mutations. In your README.md for this assignment, answer the following questions:
#What is the “size” of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
#Is this relationship significant? How do you know?

paternal_model=smf.ols(formula="paternal_dnm~1+Father_age", data=concat_DF)
#running the ordiinary least squares and storing as a variable
paternal_results=paternal_model.fit()
print(paternal_results.summary())
#print("paternal p-value is:", paternal_results.pvalues[1])








    #STEP 2.4
#Using your results from step 2.3, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband’s time of birth. 
#Record your answer and your work (i.e. how you got to that answer) in your README.md.
#predict paternal dnms for man who is 50.5 years old 


paternal_DF_dnm=concat_DF.loc[:,"paternal_dnm"]
#print(paternal_DF)
paternal_DF_age=concat_DF.loc[:,"Father_age"]
#print(paternal_DF)
concat_paternal=pd.concat([paternal_DF_age, paternal_DF_dnm], axis=1, join="inner")
#print(concat_paternal)


paternal_alone_model=smf.ols(formula="paternal_dnm~1+Father_age", data=concat_DF).fit()
predict_paternal_DNMs = pd.DataFrame({"Father_age": [50.5]})
#print(predict_paternal_DNMs)


print("the predicted number of paternal DNMs for a proband with a father who was 50.5 years old at proband TOB is:", paternal_alone_model.predict(predict_paternal_DNMs))






    #STEP2.5

#Next, you’re curious whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. 
#Using matplotlib, plot the distribution of maternal DNMs per proband (as a histogram). 
#In the same panel (i.e. the same axes) plot the distribution of paternal DNMs per proband. 
#Make sure to make the histograms semi-transparent so you can see both distributions. Upload as ex2_c.png in your submission directory.
#print(concat_DF)


from pylab import xticks

fig3, ax=plt.subplots()

number_of_maternal_DNMs=concat_DF.loc[:,"maternal_dnm"]
number_of_paternal_DNMs=concat_DF.loc[:,"paternal_dnm"]
ax.set_title( "Distributions of maternally and paternally inherited DNMs" )
ax.set_xlabel("Number of inherited DNMs per proband")
ax.set_ylabel("Frequency")
ax.hist(number_of_paternal_DNMs, alpha=0.3, label="Paternal")
ax.hist(number_of_maternal_DNMs,alpha=0.4, label="Maternal")
ax.legend()


plt.show()




    #STEP 2.6
#Now that you’ve visualized this relationship, you want to test whether there is a significant difference between the number of maternally vs. 
#paternally inherited DNMs per proband. What would be an appropriate statistical test to test this relationship? 
#Choose a statistical test, and find a Python package that lets you perform this test. If you’re not sure where to look, 
#the stats module from scipy (more here) provides tools to perform several different useful statistical tests. After performing your test, 
#answer the following answers in your README.md for this assignment:

#What statistical test did you choose? Why?
#Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.
print(sps.ttest_ind(number_of_maternal_DNMs, number_of_paternal_DNMs))









fig1.savefig("ex2_a.png")
fig2.savefig("ex2_b.png")
fig3.savefig("ex2_c.png")

plt.close(fig1)
plt.close(fig2)
plt.close(fig3)


