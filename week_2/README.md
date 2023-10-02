WEEK 2 ASSIGNMENT 

Question 1.1:
The number of reads required is: 30000

I got this result by running this code in Python: 
num_reads=int(coverage*genome_len/read_len)
print("number of reads required is:", (num_reads))

Where:
coverage= 3
genome_len=1,000,000
read_len= 100

You could also have run this code in R to find the same answer:
start_position=np.random.randint(low=0, high=high+1, size=num_reads) 
print("number start position is:", len(start_position))
number start position is: 30000

==========================================================================================================================

Question 1.3:
Using your results from Step 1.2, answer the following questions in your README.md:

In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
4.6868% of the genome has not been sequenced. 

How well does this match Poisson expectations? How well does the normal distribution fit the data?
The simulation matches the Poisson expectations very well. The normal distribution does not fit the data well. 

==========================================================================================================================

Question 1.4:
In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
0.0052% of the genome has not been sequenced. 

How well does this match Poisson expectations? How well does the normal distribution fit the data?
The simulation matches the Poisson expectation very well, similarly to part 1.3. The normal distribution still does not fit the data that well (especially in comparison to the Poisson), but the normal distribution fits the data much better here than it did in part 1.3. 

==========================================================================================================================

Question 1.5:
In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
0.0005% of the genome has not been sequenced. 

How well does this match Poisson expectations? How well does the normal distribution fit the data?
The simulation matches the Poisson expectation very well, similarly to both 1.3 and 1.4. The normal distribution fits the simulation very well, MUCH better than either 1.3 or 1.4 (the Poisson is still slightly better, though). 

==========================================================================================================================


