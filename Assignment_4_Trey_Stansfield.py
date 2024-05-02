#imports packages
import numpy as np
import pandas as pd
import math

#this reads the count txt as a dataframe
counts = pd.read_csv('test_files/argR-counts-matrix.txt', header=None,sep='\t')

#this uses the psuedo count method to add one to each cell in the matrix 
#to remove 0s
for x in range(len(counts)):
    for j in range(2,20):
        counts.loc[x,j] = counts.loc[x,j]+1

print('Pusedocount Matrix')
print(counts)

#creates a copy of the counts matrix to be used a frequency matrix
frequency = counts.copy(deep=True)

#type casts every column as a float to avoid warnings
for j in range(2,20):
    frequency[j] = frequency[j].astype('float64')

#this turns each cell into its frequency. The count/31 which is the number of 
#sequences + 4 to account for puesdo count
for x in range(len(frequency)):
    for j in range(2,20):
        frequency.loc[x,j] = frequency.loc[x,j]/31
        
print('Frequency Matrix')
print(frequency)

#creates a copy of the frequency matrix to be used a log odds matrix        
log_odds = frequency.copy(deep=True)

#This turns each frequency into its log odds. Being the frequency divided by
#the background frequency the log2 is taken of the frequency.
for x in range(len(log_odds)):
    for j in range(2,20):
        log_odds.loc[x,j] = math.log2(log_odds.loc[x,j]/0.25)
        
print('Log-Odds Matrix')
print(log_odds)

log_odds.to_csv('pssm.csv',header=None)

#empty lists of the gene IDs and score for the gene ID
genes = []
scores =[]
seqs = []

#cnt varible corresponding to each line
cnt = 0

#reads through the upstream gene sequences
with open('test_files/E_coli_K12_MG1655.400_50','r') as read:
    for x in read:
        #this prints the current line number
        #This program takes around 5 minutes to run so this lets the user
        #Know the program is running
        print(f"Line number: {cnt}")
        #This extracts the gene ID and sequence from each gene
        gene = x.split('\\')[0].strip(' ')
        seq = x.split('\\')[1].strip(' ')
        #This resets the max score and max gene per each line
        max_score = 0
        max_gene = ''
        max_seq = ''
        #this takes the sequence from each gene and reads it in an 18 bp
        #sliding window.
        for i in range(len(seq)-17):
            #this resets the local score for each sliding window
            score = 0
            #Sets the motif for the current window
            motif = seq[i:i+18]
            #reads the window bp by bp
            for j in range(18):
                #uses the log odds matrix to score based on what each bp is
                #the local score is updated iteratively
                if motif[j] == 'a':
                    score += log_odds.loc[0,j+2]
                elif motif[j] == 'c':
                    score += log_odds.loc[1,j+2]
                elif motif[j] == 'g':
                    score += log_odds.loc[2,j+2]
                elif motif[j] == 't':
                    score += log_odds.loc[3,j+2]
            #if the local window score is larger then the max for the current
            #gene then it replaces the current max score with the new one and
            #stores the associated gene ID as max gene
            if score > max_score:
                max_seq = motif
                max_gene=gene
                max_score=score
        #the max score and its associated gene ID for the gene are stored in 
        #is a list.
        seqs.append(max_seq)
        genes.append(max_gene)
        scores.append(max_score)
        #interates the cnt
        cnt += 1    

#creates a data frame of the max scores for each gene.
d = {'Gene': genes, 'Sequence': seqs,'Score': scores}
all_scores = pd.DataFrame(data=d)

#this extracts the top 30 highest storing genes
top_scores = all_scores.nlargest(30, 'Score',keep="all")

#this resets the index
top_scores = top_scores.reset_index(drop=True)

#this writes the top 30 genes to a csv
top_scores.to_csv('top_scoring_geneIDs.csv',sep=',',header=True)
