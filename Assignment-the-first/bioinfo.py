#!/usr/bin/env python

# Author: Surbhi Nahata surbhin@uoregon.edu

# Check out some Python module resources:
# - https://docs.python.org/3/tutorial/modules.html
# - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
# - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.1" # Read way more about versioning here:
# https://en.wikipedia.org/wiki/Software_versioning
DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')
Fasta_file = str(f'>header1\nAAAATTTGG\nCCCAA\n>header2\nTTTTTCCCC')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score''' 
    ord_letter = ord(letter)
    phred_letter = ord_letter - 33
    return phred_letter

def qual_score(phred_score: str) -> float:
    '''Calculates average quality score of the whole phred string'''
    cumulative_phredscore = 0
    for encode in phred_score:
        cumulative_phredscore += convert_phred(encode)
    average_score = cumulative_phredscore/len(phred_score)
    return average_score

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case
    insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA_bases):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA_bases), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA_bases = DNA_bases.upper()
    return (DNA_bases.count("G")+DNA_bases.count("C"))/len(DNA_bases)

# def oneline_fasta():
#     '''Converts multiple sequence line into one'''
#     sequence: str = "" #an empty string
#     header: str = "" #an empty header
#     i = 0
#     with open (inputfile, 'r') as fh, open (outputfile, 'w') as wr:
#         for line in fh:
#             i+=1
#             line=line.strip('\n')
#             #print(line)
#             if line.startswith(">") and header == '': 
#                 header = line
#                 #print(header)
#                 wr.write(header+"\n")
#             elif line.startswith(">"):
#                 header = line
#                 wr.write(sequence+"\n")
#                 wr.write(header+"\n")
#                 sequence = ''
#             else: 
#                 sequence+=line
#         wr.write(sequence)
# oneline_fasta()
    

def calc_median(my_list: list) -> int: 
    '''returns the median in one directional list'''
    median: int
    lstLen = len(my_list) # for the length of my_list
    index = (lstLen - 1) // 2 # index should be discrete and not decimal that's why floor division
    # imagine 9/2 = 4.5 that can't be a value for median hence 9//2 gives us 4
    if (lstLen % 2 == 0): # to check even or odd 
        median = (my_list[index] + my_list[index + 1])/2.0 # even, get the middle two values
    else:
        median = my_list[index] # odd, just get the middle value
    return median

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for
    convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    qual_score
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    print("You calcluated the correct average phred score")
    
    validate_base_seq
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

    gc_content
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("correctly calculated GC content")

    calc_median
    assert calc_median([1,2,3]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    print("Median successfully calculated")

    # oneline_fasta
    # assert oneline_fasta("AAA", "\n", "CCC") == True, "You have not created one line seq"
    # print("You have succesfully created one line sequence")
    