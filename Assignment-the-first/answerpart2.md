1. Multiplexing is a process by which barcode or indices both in front(index1) and back(index2) of biological reads are added. We know these index location and sequence, and in an ideal world these indexes are complementary to each other and don't have any mismatches, low quality or unkowns. However, some times these indexes have N(AANTTTC), mismatch of base pairs giving us low quality(AAAAA-ABABAB), and show index hopping (where the index is of correct sequence and is present in the index library, however is incorrectly paried to each other for example index1(AAAAA)-index2(BBBBB)). The algoritham that we will devolop in this assignment will seperate files based on finding the mismatch, unknown, index hopping and matched index read. As we only want files that have matched index read. The other files are undesirable or we can correct those down the line but for this assignment we won't be doing that. 
2. Output is going to have:
a) Unknown/mismatch/N's/low quality (2 files- read 1 and read 2): A file where the index is unknown and has indexes that are mismatched and have N's in the sequence.
b) Index hopping (2 files -read 1 and read2): If there is index hopping i.e indexes are known but are mismatch to each other
c) Matched files (48 files - read 1 and read2): These files have matched indexes i.e both index 1 and index 2 are known and reverse complementary to each other.
3.  DNA = 'ACGTCC'
def reverse_complementary(str) ->str:

    '''Takes a string of nucleotide and swaps A for T, T for A, G for C and C for G by using numpy to reverse the string'''
    
    Reverse_DNA = reverse complements ACGTCC to GGACGT 
    return Reverse_DNA
