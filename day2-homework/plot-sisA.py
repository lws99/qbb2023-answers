#!/usr/bin/env python

import matplotlib.pyplot as plt 
import numpy as np



transcripts = np.loadtxt( "all_annotated.csv", delimiter=",", usecols=0, dtype="<U30", skiprows=1 )
#print( "transcripts: ", transcripts[0:5] )

samples = np.loadtxt( "all_annotated.csv", delimiter=",", max_rows=1, dtype="<U30" )[2:]
#print( "samples: ", samples[0:] )

data = np.loadtxt( "all_annotated.csv", delimiter=",", dtype=np.float32, skiprows=1, usecols=range(2, len(samples) + 2) )
#print( "data: ", data[0:5, 0:5] )


# Find row with transcript of interest
for i in range(len(transcripts)):
    if transcripts[i] == 'FBtr0073461': 
        row = i
        #print(row)



# Find columns with samples of interest
cols_female = []
cols_male=[]
for i in range(len(samples)):
    if "female" in samples[i]:
        cols_female.append(i)
    else:
    	cols_male.append(i)


#print(cols_female)
#print(cols_male)


# Subset data of interest
expression_female = data[row, cols_female]
expression_male=data[row, cols_male]

#print(expression_female)
#print(expression_male)


# Prepare data
x_female = samples[cols_female]
y_female= expression_female

print(x_female)
#print(expression_female)

x_male=samples[cols_male] #we technically do not need this because the female and the males will use the same x axis for plotting
y_male=expression_male 
print(x_male)

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


ax.set_xticklabels(labels=labels, rotation=90) #rotates the x axis 90 degrees and labels with dev stage
ax.legend()

plt.show()
fig.savefig("sisA FBtr0073461.pdf" )
plt.close( fig )
