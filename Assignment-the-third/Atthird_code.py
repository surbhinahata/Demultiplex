#!/usr/bin/env python
# author: Surbhi Nahata
# Date: 8/2/2023

###########
# imports #
###########

#import bioinfo #for convert phred function
import argparse #argparse to make the code general for any file
import gzip #Use gzip to access file content without decompressing it (decompressing takes up a lot of memory so avoid it)
import matplotlib.pyplot as plt #to build graphs
import math 
#import itertools
import numpy

#############
# Arguments #
#############

def get_args():
    parser = argparse.ArgumentParser(description="Gets fastq file for index demultiplexing")
    # You want to take for input files and read it at the same time read1,read2,read3 and read4
    parser.add_argument("-bf1", "--read1", help = "biological read 1 with fastq data", required = True)
    parser.add_argument("-bf2", "--read2", help = "biological read 2 with fastq data", required = True)
    parser.add_argument("-if1", "--index1", help = "index 1 with fastq data", required = True)
    parser.add_argument("-if2", "--index2", help = "index 2 with fastq data", required = True)
    parser.add_argument("-u1", "--unknown1_file", help = "output file name", required = True)
    parser.add_argument("-u2", "--unknown2_file", help = "output file name", required = True)
    parser.add_argument("-h1", "--hopped1_file", help = "output file name", required = True)
    parser.add_argument("-h2", "--hopped2_file", help = "output file name", required = True)
    parser.add_argument("-c", "--cutoff", help = "Takes a Qscore cut off value", type= int, required = True)
    args = parser.parse_args()
    return args

args = get_args()
Read1 = args.read1
Read2 = args.index1
Read3 = args.index2
Read4 = args.read2
cutoff = args.cutoff
unknownoutR1 = args.unknown1_file
unknownoutR2 = args.unknown2_file
indexhoppingR1 = args.hopped1_file
indexhoppingR2 = args.hopped2_file

###########
# Globals #
###########

thefile = "/projects/bgmp/surbhin/bioinfo/Bi622/Demultiplex/Assignment-the-third/indexes.txt" #fix it to argsparse if you get time since we want to keep it general and not hardcode it 

############# 
# Functions #
#############

def reverse_complement(str):
    '''A function that takes ATGC nucleotide and swaps it to GCAT with reverse complementation '''
    
    comp_dict: dict = {"A":"T", "T":"A","C":"G","G":"C","N":"N"} #creating dictionary that holds complementary nucleotide
    compstring: str = "" #empty string to which you're going to add complementary nucleotide after looping over nucleotide in DNA_input

    for nt in seq3: #for nucleotide in DNA_input sequence
        compstring += comp_dict[nt] #keep adding the dictinary value into the compstring so TGCAA...
    return compstring[::-1]

def getfastqfile (filehandle) -> tuple[str, str, str, str]:

    '''A function that takes multiple file at the same time and splits them into header, sequence, comment and qscore'''

    header = filehandle.readline().strip() # gets header
    seq = filehandle.readline().strip() # gets sequence
    comment = filehandle.readline().strip() # gets comment
    qscore = filehandle.readline().strip() # gets quality score
    return (header, seq, comment, qscore)

def convertindexset (thefile: str) -> set[str]: #function can only take varible name in the paranthesis you cant give it a file name that's why defined it earlier
    '''A function that gets the file containing index and spits out the index sequence in return of function'''
    index_set = set() #an empty set
    with open(thefile, 'r') as fh:
        fh.readline()
        for line in fh:
            line = line.strip("\n").split("\t")
            line = line[-1] #taking out the 5th element in the indexes file which is the index sequence
            #print(line)
            index_set.add(line) #method which add stuff to the set
    return index_set

def openindexfile(index_set: set):
    '''Given a list of string indexes, open up a couple fast q of files for each index, 
    storing those file handles into a dictionary that I return.'''

    fh_dict: dict = {} #an empty dictionary that will store index/barcode as the key and value as filehandle1 and filehandle 2

    for seq in index_set: #it will loop over the seq in index
        fh1 = open (f"{seq}_R1.fastq", 'w')
        fh2 = open (f"{seq}_R2.fastq", 'w')
        fh_dict[seq] = (fh1, fh2) #tuple because we know what the value is and won't change
    return fh_dict

########
# Main #
########

index_set = convertindexset(thefile)
fh_dict = openindexfile(index_set)

#opening the write up file above the loop is very important. If you write it under the loop it will take millions of years to open and then write stuff
with open(unknownoutR1, 'w') as wr1, open(unknownoutR2, 'w') as wr2, open(indexhoppingR1, 'wt') as hwr1, open(indexhoppingR2, 'wt') as hwr2: 
    with gzip.open(Read1, 'rt') as R1, gzip.open(Read2, 'rt') as I1, gzip.open(Read3, 'rt') as I2, gzip.open(Read4, 'rt') as R2:
        unknown_counter = 0
        Total_counter = 0
        index_hopped = {} #empty dictionary to count the number of index hopped combos
        matched_count = {} #empty dictionary to count the matched numbers 
        
        while True:
            header1, seq1, comment1, qscore1 = getfastqfile(R1) #in the read file get header 1, sequence 1, +, and qscore1 for biological read 1 
            header2, seq2, comment2, qscore2 = getfastqfile(I1) #same as above but for index1
            header3, seq3, comment3, qscore3 = getfastqfile(I2) #same as above but for index2
            header4, seq4, comment4, qscore4 = getfastqfile(R2) #same as aboe but for biological read 4
            if header1 == "": #break if you find empty line
                break
            Total_counter+=1
            # quick test to see if getfastqfile is working
            #wr.write(f'{header1}\n{seq1}\n{comment1}\n{qscore1}\n')
            #wr.write(f'{header2}\n{seq2}\n{comment2}\n{qscore2}\n')
            #wr.write(f'{header3}\n{seq3}\n{comment3}\n{qscore3}\n')
            #wr.write(f'{header4}\n{seq4}\n{comment4}\n{qscore4}\n')

            rev_seq3 = reverse_complement(seq3)
            #print(reversed)
        
            header1 += " " + seq2 + '-' + rev_seq3
            header4 += " " + seq2 + '-' + rev_seq3
            
            if ("N" in seq2) or ("N" in rev_seq3): #if N present in the index
                wr1.write(f'{header1}\n{seq1}\n{comment1}\n{qscore1}\n') # read 1
                wr2.write(f'{header4}\n{seq4}\n{comment4}\n{qscore4}\n') # read 1
                unknown_counter+=1
            
            elif seq2 not in fh_dict or rev_seq3 not in fh_dict: #if index not in dictionary 
                wr1.write(f'{header1}\n{seq1}\n{comment1}\n{qscore1}\n') # read 1
                wr2.write(f'{header4}\n{seq4}\n{comment4}\n{qscore4}\n') # read 1
                unknown_counter+=1
                
            elif seq2 == rev_seq3:
                fh_dict[seq2][0].write(f'{header1}\n{seq1}\n{comment1}\n{qscore1}\n') #I opened the fh_dict in the function already. key seq 2 is getting the 0th value (see function for the syntax)
                fh_dict[seq2][1].write(f'{header4}\n{seq4}\n{comment4}\n{qscore4}\n') #this is getting the key seq 2 and getting the 1st value
                if seq2 in matched_count:
                    matched_count[seq2]+=1
                else:
                    matched_count[seq2]=1

            else:
                hwr1.write(f'{header1}\n{seq1}\n{comment1}\n{qscore1}\n')
                hwr2.write(f'{header4}\n{seq4}\n{comment4}\n{qscore4}\n')
                if (seq2, rev_seq3) in index_hopped: #if in dictionary add into it 
                    index_hopped[(seq2, rev_seq3)]+=1
                else:
                    index_hopped[(seq2, rev_seq3)]=1

print(Total_counter)
# print(unknown_counter)
# print((unknown_counter/Total_counter)*100)
# print(index_hopped)
print(matched_count)

#Total_count=sum(index_hopped.values())+sum(matched_count.values())+unknown_counter
#print(Total_count)

with open('indexhopped_data.tsv','w') as w1, open('matchedindex_data.tsv','w') as w2, open('complete_data.tsv','w') as w3:
    print(f'index-pair\thopped count\tpercent of total hopped pairs',file=w1)
    print(f'index\tmatched count\tpercent index in sample\tpercent index overall',file=w2) 

    for indexed_pair in index_hopped: #writing the index pair to a file that will count each instance
        percent_hopped=index_hopped[indexed_pair]/Total_counter*100
        print(f'{indexed_pair}\t{index_hopped[indexed_pair]}\t{percent_hopped}%',file=w1)

           
    for indexm in matched_count: #writing how many indexes were matched
        percent_matched=matched_count[indexm]/sum(matched_count.values())*100 #calculating the number of unique indexes in matched 
        percent_matched_total=matched_count[indexm]/Total_counter*100 #calculating the number of unique matched indexes in total
        print(f'{indexm}\t{matched_count[indexm]}\t{percent_matched}%\t{percent_matched_total}%',file=w2)
            
    print(f'type of data\tcount\tpercent of total',file=w3)
    percent_hopped=sum(index_hopped.values())/Total_counter*100
    percent_matched=sum(matched_count.values())/Total_counter*100
    percent_unknown=unknown_counter/Total_counter*100
    print(f'hopped\t{sum(index_hopped.values())}\t{percent_hopped}%', file=w3)
    print(f'matched\t{sum(matched_count.values())}\t{percent_matched}%', file=w3)
    print(f'unknown\t{unknown_counter}\t{percent_unknown}%', file=w3)
