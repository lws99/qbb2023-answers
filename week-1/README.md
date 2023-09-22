# Homework N Answers:
## Question 2.2:
What is the “size” of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 2.1?
As maternal age increases by one year, so does the number of maternal DNMs by 0.3776. This increase in maternal DNMs is the size of this relationship between maternal age and maternal DNMs. This is similar to the plot I got in 2.1 where there was a slight positive correlation between maternal age and number of maternal DNMs. 


Is this relationship significant? How do you know?
Yes the relationship is significant because the p-value of the relationship is 6.878208052158821e-24.


## Question 2.3:
What is the “size” of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
Is this relationship significant? How do you know?
As paternal age increases by one year, so does the number of paternal DNMs by 1.3538. This increase in paternal DNMs is the size of this relationship between paternal age and paternal DNMs. This is similar to the plot I got in 2.1 where there was a clear positive correlation between paternal age and number of paternal DNMs. This relationship was also more obvious than the relationship that I saw between maternal age and maternal DNMs, which is reflected in the paternal and maternal effect sizes. 


Is this relationship significant? How do you know?
Yes the relationship is significant because the p-value of the relationship is 1.552293615688712e-84.



## Question 2.4:
Using your results from step 2.3, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband’s time of birth. Record your answer and your work (i.e. how you got to that answer) in your README.md.

My model predicted that a proband with a father who is 50.5 years old at the time of proband's birth would have 78.695457 paternal DNMs. 

This is how I got my answer:
1. Create a new dataframe that only contains father DNMs and father age as columns
paternal_DF_dnm=concat_DF.loc[:,"paternal_dnm"]
paternal_DF_age=concat_DF.loc[:,"Father_age"]
concat_paternal=pd.concat([paternal_DF_age, paternal_DF_dnm], axis=1, join="inner")

2. Create a new dictionary that contains the paternal age that I want to predict paternal proband number from. 
predict_paternal_DNMs = pd.DataFrame({"Father_age": [50.5]})

3. Use the linear regression model that I created that predicts paternal DNM number from father's age to find the predicted number of paternal probands for a father who is 50.5.
paternal_alone_model=smf.ols(formula="paternal_dnm~1+Father_age", data=concat_DF).fit()
print("the predicted number of paternal DNMs for a proband with a father who was 50.5 years old at proband TOB is:", paternal_alone_model.predict(predict_paternal_DNMs))



##  Question 2.6:
What statistical test did you choose? Why?
Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.

I chose to use an unpaired t-test because the two groups that the DNMs come from (mother and father) are independent from each other. 
The t-test that I ran was siginificant. Thus, there are more paternally inherited DNMs than maternal DNMs. 
