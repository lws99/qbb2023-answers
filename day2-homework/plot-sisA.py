#!/usr/bin/env python

import matplotlib.pyplot as plt 
import numpy as np



transcripts = np.loadtxt( "all_annotated.csv", delimiter=",", usecols=0, dtype="<U30", skiprows=1 )
#print( "transcripts: ", transcripts[0:5] )


#EXERCISE 2
dataset_open=open("all_annotated.csv")
dataset=dataset_open.readlines()
#print(dataset)


				#only need to read out the transcripts 

transcripts_loop=[]

for file in dataset:
	file = file.rstrip()
	file_list=file.split(',')
	transcripts_loop.append(file_list[0])

transcripts_final = transcripts_loop[1:]
#print(transcripts_final)


samples = np.loadtxt( "all_annotated.csv", delimiter=",", max_rows=1, dtype="<U30" )[2:]
#print( "samples: ", samples[0:] )

data = np.loadtxt( "all_annotated.csv", delimiter=",", dtype=np.float32, skiprows=1, usecols=range(2, len(samples) + 2) )
#print( "data: ", data[0:5, 0:5] )




#EXERCISE 1

# Find row with transcript of interest
for i in range(len(transcripts)):
    if transcripts[i] == 'FBtr0073461': 
        row = i
        #print(row)





# Find columns with samples of interest and make 2 lists, one containing the male data and containing the female data
cols_female = []
cols_male=[]
for i in range(len(samples)):
    if "female" in samples[i]:
        cols_female.append(i)
    else:
    	cols_male.append(i)


#print(cols_female)
#print(cols_male)






# Subset data of interest of males and females
expression_female = data[row, cols_female]
expression_male=data[row, cols_male]

#print(expression_female)
#print(expression_male)






# Prepare female data
x_female = samples[cols_female]
y_female= expression_female

#print(x_female)
#print(y_female)
#print(expression_female)

#Prepare male data
x_male=samples[cols_male] #we technically do not need this because the female and the males will use the same x axis for plotting
y_male=expression_male 

#print(x_male)
#print(y_male)

#Prepare 2* male data
twice_male_y_data=2*(np.array(y_male))
#print(twice_male_y_data)



# Plot data
plt.rcParams['text.usetex'] = False # gives the ability to italicize
labels=['10', '11', '12', '13', '14A', '14B', '14C', '14D']
fig, ax = plt.subplots()
ax.set_title("$\it{sisA}$" ) #italicizing the title
ax.set_xlabel("developmental stage")
ax.set_ylabel("mRNA abundance (RPKM)")
ax.set_ylim(0,250)

ax.plot(x_female, y_female, c = "red", label="Female")
ax.plot(x_female, y_male, c = "blue", label="Male")
ax.plot(x_female, twice_male_y_data, c = "lightblue", label="2*Male")


ax.set_xticklabels(labels=labels, rotation=90) #rotates the x axis 90 degrees and labels with dev stage
ax.legend()

plt.show()
fig.savefig("sisA FBtr0073461.pdf" )
plt.close( fig )


