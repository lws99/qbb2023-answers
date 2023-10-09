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

Question 2.4:
Now, use dot to produce a directed graph. Record the command you used in your READMD.md.


cat graph.txt | dot -Tpng > ex2_digraph.png

============================================================================================================================

Question 2.5:

Assume that the maximum number of occurrences of any 3-mer in the actual genome is five. Using your graph from Step 2.4, write one possible genome sequence that would produce these reads. Record your answer in your README.md.



TTCATTCTTATTGATTT


=============================================================================================================================

Question 2.6:

In a few sentences, what would it take to accurately reconstruct the sequence of the genome? 

You would need to know where the start of the genome is in the directed graph. You would also need to know how many occurances (repeats) there are for EACH kmer in the genome. Understanding the repeats allows you to understand what direction to take in the directed graph when you have multiple directional choices to the next kmer. This directionality, as determined by repeats for EACH kmer, matters because it determines the final sequence that you end up with.


=====================================================================================================================================

Question 3.1:

Use ChatGPT (or Bard or your favorite LLM) to write an essay on why you are interested in genomics. Make sure to ask for references. 





Q: Write an essay on why I should be interested in genomics. Cite your sources. 




A:
Introduction

Genomics, the study of an organism's complete set of genes and their functions, has emerged as a groundbreaking field with profound implications for individuals and society as a whole. This essay aims to elucidate why you should be interested in genomics by highlighting its multifaceted benefits and its capacity to reshape healthcare, advance scientific knowledge, and address ethical and societal concerns. Citing credible sources, we will explore the compelling reasons to engage with this transformative discipline.

Personalized Medicine and Health Benefits
One of the most compelling reasons to be interested in genomics is its pivotal role in advancing personalized medicine. Genomic information allows healthcare providers to tailor treatments to an individual's genetic makeup, ensuring more effective and safer interventions. As reported in a study published in Nature Reviews Genetics (Pirmohamed et al., 2013), pharmacogenomics, a branch of genomics, has the potential to revolutionize healthcare by optimizing drug prescriptions and minimizing adverse reactions.

Disease Prevention and Early Detection
Genomics empowers individuals to take proactive measures in preventing genetic diseases. By understanding their genetic predispositions, individuals can make informed lifestyle choices and undergo early screenings, as noted in the New England Journal of Medicine (Khera et al., 2018). This knowledge can significantly reduce the risk of developing certain conditions, leading to healthier lives.

Revolutionizing Cancer Research
Cancer research has been transformed by genomics, which enables the identification of specific genetic mutations driving cancer and the development of targeted therapies. A seminal article in Nature (Garraway & Lander, 2013) underscores how genomics has reshaped our approach to understanding and treating cancer, potentially leading to more effective treatments and improved survival rates.

Ethical and Societal Implications
Engaging with genomics is essential to navigate its ethical and societal implications. The field raises questions about genetic privacy, discrimination, and the moral quandaries surrounding gene editing technologies like CRISPR-Cas9. The report by the National Academies of Sciences, Engineering, and Medicine examines these critical ethical considerations, emphasizing the importance of informed public discourse.

Applications Beyond Medicine
Genomics extends its influence beyond healthcare. In agriculture, for instance, it plays a pivotal role in developing crops with improved traits such as disease resistance and increased nutritional value, as acknowledged by the Food and Agriculture Organization of the United Nations. This demonstrates how genomics can address global challenges like food security.

Personal Ancestry and Genealogy
For those with a curiosity about their origins and heritage, genomics offers the opportunity to trace ancestry and genealogy through commercial genetic testing services like 23andMe and AncestryDNA. These services provide individuals with insights into their genetic heritage and connections with relatives they may not have known about.

Scientific Advancement and Intellectual Stimulation
Lastly, for those with a broader interest in science and technology, genomics provides exciting opportunities for research and discovery. This multidisciplinary field combines biology, genetics, computer science, and more, offering a stimulating avenue for exploration and innovation.

Conclusion

In conclusion, genomics is not merely a scientific endeavor; it is a transformative force that impacts every facet of our lives. It has the potential to enhance our health, reshape medical treatments, fuel scientific discovery, and provoke critical ethical discussions. By staying informed and engaged with genomics, you can be part of a revolution that holds the promise of a healthier, more informed, and ethically conscious future.

As with any rapidly evolving field, it is essential to consult reputable sources, scientific journals, government agencies, and academic institutions to stay up-to-date with the latest developments in genomics. Embracing genomics means embracing the future of healthcare, science, and society as a whole.




=====================================================================================================================================

Step 3.2:

In your README.md, comment on the output from the LLM: Does it make logical sense? Does it include any phrases you would not have written? Do the cited papers exist and support the claims from the LLM?



The output is mostly logical, but fairly redundant. For example, Chat GPT goes into both personalized medicine and personalized cancer treatments in detail. This could have been one point where cancer treatment personalization is an example of personalized medicine. The response by Chat GPT is generally very vauge, especially in the point about the ethical and societal implications. Chat GPT just broadly references some report by the National Academies of Sciences, Engineering, and Medicine that I am unfamiliar with and cannot find because there is no direct source cited (same thing in the applications beyond medicine). I would have cited a specific source here and not brought up the information as if the reader already knew it. In the conclusion, I do like the fact that Chat GPT provides general places for the user to continue looking for information on genomics, but I wish that Chat GPT had provided some specific examples of reputable places to continue looking into genomics. 

Overall, I would have been a lot more specific with my essay than Chat GPT was, but there are pieces of good information in what Chat GPT provided. I also don't think that the ethical and societal implications paragraph fit well into the essay that Chat GPT wrote. It is an important point, but it feels out of place because of the lack of context that Chat GPT provided.

Citation accuracy:

The citation New England Journal of Medicine (Khera et al., 2018) has been cited incorrectly, as I am unable find a paper from Khera et al 2018 in the NEJM. I did find one paper from Khera et al. related to genetic risk of coronary heart disease from 2016 in the NEJM.

The citation Nature Reviews Genetics (Pirmohamed et al., 2013) is also incorrect. There is a scientist by this name, but they have never published anything as first author in Nature Reviews Genetics, and they did not publish anything in this journal in 2013. 


All other citations are extremely vague and Chat GPT essentially assumes that the reader already knows the article it is citing, which is a bad assumption to make. 











