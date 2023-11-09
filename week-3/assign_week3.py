#!/usr/bin/env python

import numpy as np
import sys


import pandas as pd
from fasta import readFASTA

"""
Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:

A FASTA-style file containing two sequences to align
A text file containing the scoring matrix you’d like to use for this alignment
The penalty for gaps in your alignment (so if users wanted to penalize gaps by subtracting 10 from the alignment score, they would input -10)
The filepath to write your alignment to
Additionally, your script should print out the number of gaps in the first sequence, the number of gaps in the second sequence, and the score of the final alignment.

You’ll run your script twice:

Align the CTCF DNA transcript sequences from human and mouse using the HOXD70 scoring matrix and a gap penalty of -300.
Align the CTCF amino acid sequences from human and mouse using the BLOSUM62 scoring matrix and a gap penalty of -10.
NOTE: The DNA sequences are fairly long, and as such the DNA alignment may take a few minutes to run. 
#We recommend testing your code with the protein alignment first (or even just a couple of small test sequences), 
#and then running the DNA alignment when you’re confident it’s working.


"""

#============================================================================================================
#                #Step 1.1: Read in your parameters

#Use sys.argv to read in the parameters you need to run your script from the command line. When reading in the fasta file, 
#you can use readFASTA to store the sequence ids and sequences as follows:

#For the scoring matrix, it probably makes sense to store it in a pandas dataframe.

#HINT: Uh oh! Looks like the scoring matrices don’t use a consistent field separator… Maybe pd.read_csv() 
#has an argument that lets you separate columns based on an arbitrary amount of whitespace…



#put into the command line: python assign_week3.py CTCF_38_M27_AA.faa <-- (filename)
FASTA=sys.argv[1]
scoring_matrix_sys=sys.argv[2]
#sys.argv always reads in data as a str so you need to convert it to a float
gap_penalty=float(sys.argv[3])
output_file=sys.argv[4]


#FOR RUNNING THE AA DATASET
#./assign_week3.py CTCF_38_M27_AA.faa  BLOSUM62.txt -10 x



input_sequences = readFASTA(open(FASTA))


seq1_id, sequence1 = input_sequences[0]


seq2_id, sequence2 = input_sequences[1]


scoring_matrix= pd.read_csv(scoring_matrix_sys, delim_whitespace= True)
#scoring_matrix=float(scoring_matrix)



#===========================================================================================================================================

            #Step 1.2: Initializing matrices

#You’ll need two matrices to carry out the Needleman-Wunsch algorithm: an F-matrix that stores the score of each “optimal” sub-alignment 
#(this is the one we created in class), as well as a traceback matrix that allows you to determine the optimal global alignment (as a path through this matrix). 
#Initialize two empty matrices for these purposes.

#HINT: With sequence 1 of length m and sequence 2 of length n, both matrices should be of size (m+1)×(n+1), to account for potential leading gaps in either sequence.

#np.zeros need to be fed in as strings for the tb matrix arg is --> ,str

#Initialize two empty matricies

#final matricies will need to be the same size the +1 accounts for if the first two bases do not match on each sequence (this is not always the case) --> -1 -1 in d of the for loop for f matrix undoes this +1 +1 if you don't have a leading gap 
F_matrix=np.zeros((len(sequence1)+1, len(sequence2)+1))

tb_matrix=np.zeros((len(sequence1)+1, len(sequence2)+1), str)




#manually fill in first two positions of the alignment and then the algorithm does the rest 

for i in range(len(sequence1)+1):
    F_matrix[i,0]=i*gap_penalty



for j in range(len(sequence2)+1):
    F_matrix[0,j]=j*gap_penalty




for i in range(len(sequence1)+1):
    tb_matrix[i,0]='v'





for j in range(len(sequence2)+1):
    tb_matrix[0,j]='h'



#==============================================================================================================================================================

                #Step 1.3: Populating the matrices

#Follow the steps of the needleman-wunsch algorithm discussed in class to populate the two matrices.

#When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, 
#resolve the tie in the order (aligning -> gap in sequence 1 -> gap in sequence 2).

# print(seq1_id, sequence1)
# print(len(sequence1))
# print(len(sequence2))
# print(F_matrix.shape)

for i in range(1, len(sequence1)+1):
    for j in range(1, len(sequence2)+1):
        #print(i, j)

        d = F_matrix[i-1, j-1] + scoring_matrix.loc[sequence1[i-1], sequence2[j-1]]
        h = F_matrix[i,j-1] + gap_penalty
        v = F_matrix[i-1,j] + gap_penalty

        F_matrix[i,j] = max(d,h,v)

        if d == max(d,h,v):

            tb_matrix[i,j] = 'd'
        elif v == max(d,h,v):

            tb_matrix[i,j] = 'v'

        elif h == max(d,h,v):

            tb_matrix[i,j] = 'h'

        else:
            print('fix your logic')




#print(tb_matrix)
#print(F_matrix)




#==================================================================================================================================
#Step 1.4: Find the optimal alignment

#Use the traceback matrix to find the optimal alignment between the two sequences. Start at the bottom right corner and follow a path 
#backwards through the traceback matrix until you reach the top left of the matrix, building the alignment as you go. You should end up with two 
#strings of the same length, one for each sequence. Gaps in the sequences should be denoted with a hyphen (-). For example, if your input sequences were 
#TACGATTA and ATTAACTTA your final alignment might look something like:

#i is row j is column

#print(sequence1)
#print(sequence2)

#each is an empty string
seq1_aligned=""
seq2_aligned=""

i = len(sequence1)
j = len(sequence2)

#print(i, j)
while i > 0 and j > 0:
    #print(tb_matrix[i,j])
    seq1_index=i-1 #each position starting at i-1
    seq2_index=j-1
    if tb_matrix[i,j] == "d":
        #seq1_aligned.append(sequence1[i])
        #seq2_aligned.append(sequence2[j])
        seq1_aligned=sequence1[seq1_index]+seq1_aligned
        seq2_aligned=sequence2[seq2_index]+seq2_aligned #order of the addition matters 
        i = i -1 #acting like a counter
        j = j -1 
    elif tb_matrix[i,j] =="h":
        # seq1_aligned.append("-") 
        # seq2_aligned.append(sequence2[j])
        seq1_aligned="-"+seq1_aligned
        seq2_aligned=sequence2[seq2_index]+seq2_aligned
        #i=i
        j=j - 1
    elif tb_matrix[i,j] == "v":
        # seq1_aligned.append(sequence1[i]) 
        # seq2_aligned.append("-")
        seq1_aligned=sequence1[seq1_index]+seq1_aligned
        seq2_aligned="-"+seq2_aligned
        i=i - 1
        #j=j


#
print(seq1_aligned)
print(seq2_aligned)
#print(sequence1)
#print(sequence2)

print("seq1 number of gaps is: ", seq1_aligned.count("-"))
print("seq2 number of gaps is: ", seq2_aligned.count("-"))



print("matrix score is: ", F_matrix[-1,-1])

output_writing=open(output_file, "w")

for position in range(0, len(seq1_aligned), 10):
    output_writing.write(seq1_aligned[position:position+10]+"\n") 
    output_writing.write(seq2_aligned[position:position+10]+"\n")

output_writing.close()


#===========================================================Step 1.5: Write the alignment to the output

# Write the alignment to the output file specified in the command line.

# ALSO, make sure your script prints out the additional requested information:

# the number of gaps in each sequence and
# the score of the alignment
# For both alignments (DNA and AA), record these values in your README.md.

















