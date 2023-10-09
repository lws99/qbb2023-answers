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
f_matrix=np.zeros((len(seq1_id)+1, len(seq2_id)+1))
tb_matrix=np.zeros((len(seq1_id)+1, len(seq2_id)+1))




#manually fill in first two positions of the alignment and then the algorithm does the rest 

for i in range(len(seq1_id)+1):
    f_matrix[i,0]=i*gap_penalty



for j in range(len(seq2_id)):
    f_matrix[0,j]=j*gap_penalty





#==============================================================================================================================================================

                #Step 1.3: Populating the matrices

#Follow the steps of the needleman-wunsch algorithm discussed in class to populate the two matrices.

#When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, 
#resolve the tie in the order (aligning -> gap in sequence 1 -> gap in sequence 2).

for i in range(1, len(seq1_id)+1):
    for j in range(1, len(seq2_id+1)):
        if sequence1[i-1] == sequence2[j-1]:  #if there is a match in the positions, then move diagonally
            d = F_matrix[i-1, j-1] + match_score
        else: 
            d = F_matrix[i-1, j-1] + mismatch_score
        h = F_matrix[i,j-1] + gap_penalty
        v = F_matrix[i-1,j] + gap_penalty

        F_matrix[i,j] = max(d,h,v)




#for a in range(len(scoring_matrix)+1):
 #   tb_matrix[i,0]=a





