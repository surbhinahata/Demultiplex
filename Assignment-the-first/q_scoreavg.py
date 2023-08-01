#!/usr/bin/env python

import bioinfo
import argparse
import gzip
import matplotlib.pyplot as plt 

def get_args():
    parser = argparse.ArgumentParser(description="Gets fastq file for finding average qscore distribution")
    parser.add_argument("-f", "--file", help = "file with fastq data", required = True)
    parser.add_argument("-l", "--seq_length", type = int, help = "Enter sequence length", required = True)
    parser.add_argument("-o", "--output_file", help = "output file name", required = True)
    args = parser.parse_args()
    return args

args = get_args()
file = args.file
seq_len = args.seq_length
out = args.output_file

def init_list(lst: list, value: float=0.0) -> list:
    
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    
    for phredvalues in range(seq_len):
        lst.append(value)
    return lst

my_list: list = []
my_list = init_list(my_list)

def populate_list(file: str) -> tuple[list, int]:

    """Stores phred score of base pair in 101 bp lane and counts the sum of each column""" # eg for test file 400xcol0, 400xcol1
    
    with gzip.open(file, 'rt') as fh:
        my_list = init_list([]) # a different way of typing my_list: list = [] my_list = init_list(my_list)
        i = 0
        for line in fh:
            i+=1
            line = line.strip('\n')
            if i%4 == 0:
                for index, encode in enumerate(line):
                    converted_score = bioinfo.convert_phred(encode)
                    my_list[index] += converted_score
    return my_list, i

my_list, num_lines = populate_list(file)

for f, value in enumerate(my_list):
    my_list[f] = (my_list[f])/(num_lines/4)

x = range(seq_len)
y = my_list

plt.bar(x, y, color='violet') 
plt.xlabel('Base index') 
plt.ylabel('Mean qscore') 
plt.title(f"mean qscore at each base for {out}") # displaying the title

#plt.show()
plt.savefig(f'{out}.png')