'''
This is going to filter out the sequence that is too short or too long (< 210 or > 410 bp)
'''
import pandas as pd
import numpy as np
import os, sys
from Bio import SeqIO

WTSEQENCE = "AGACATAAAGGGACTCACCTCATGGGTGTAAGTACGAACAGGGACTCCAAAAGCTCTGGCTAGTCCAAAGTCTGCTAGCTTGA"
NUMOFNT = len(WTSEQENCE)

path = os.getcwd()
files = os.listdir(path)
files_fasta = [f for f in files if f[-5:] == 'fasta']

for f in files_fasta:
    
    seqdic = {}
    start_dna_seq = 'GAGAAAAA'
    end_dna_seq = 'TGGCCCCC'
    startpos = 0
    endpos = 0
    seqarray = []
    numofntseq = []
    insertseq_array = []
    new_list = []
    new_list2 = []
    another_list = []
    classifier_array = []
    nt_length_array = []
    #handle fasta using biopython
    handle = open(f, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        
        seqdic[record.id] = record.seq
    handle.close()

    s1 = pd.Series(seqdic, name='seq')
    s1 = s1.to_frame(name='seq')

    for index, row in s1.iterrows():
        new_list.append(str(row[0]))

    for seq in new_list:
        another_list.append(seq.replace('-', ''))

    for clean_seq in another_list:
        nt_length_array.append(len(clean_seq))
    
    #filter out too short and too long sequence
    s1['number_nt'] = nt_length_array
    s1 = s1[s1['number_nt'] < 410]
    s1 = s1[s1['number_nt'] > 210]
    
    for index, row in s1.iterrows():
        new_list2.append(str(row[0]))
    
    for clean_seq in new_list2:
        startpos = clean_seq.find(start_dna_seq) + 8
        endpos = clean_seq.find(end_dna_seq)
        if clean_seq.find(start_dna_seq) == -1 or endpos == -1:
            numofntseq.append(0)
            insertseq_array.append("False")
        else:
            numofntseq.append(endpos - startpos)
            insertseq_array.append(clean_seq[startpos:endpos])

    s1['numofntseq'] = numofntseq
    s1['insertseq'] = insertseq_array

    for num, seq in zip(s1['numofntseq'], s1['insertseq']):
        if num == NUMOFNT and seq == WTSEQENCE:
            classifier_array.append("WT")
        elif abs(num-NUMOFNT)%3 == 0 and num !=NUMOFNT:
            classifier_array.append("in frame")
        elif abs(num-NUMOFNT)%3!= 0 and num != 0:
            classifier_array.append("out of frame")
        elif num == NUMOFNT:
            classifier_array.append("missense")
        elif seq == 'False':
            classifier_array.append("no match")
        else:
            classifier_array.append("Exception")
            
    s1['group'] = classifier_array
    #s1 = s1[s1.group != 'no match']
    s2 = s1.groupby(['group']).size().reset_index(name='count')

    total = 0.0
    for index, row in s2.iterrows():
        total += row['count']

    percent_array = []

    for index, row in s2.iterrows():
        percent_array.append('{0:.2f}%'.format(row['count']/total*100))
        
        
        

    s2['percent'] = percent_array
    
    

    writer = pd.ExcelWriter('multipleFasta{0}.xlsx'.format(f[0:3]))
    s1.to_excel(writer, 'Sheet1')
    s2.to_excel(writer, 'Sheet2')
    writer.save()
