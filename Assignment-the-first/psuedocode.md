Read1(R1) and Read2(R4) are biological reads that have 2 complementary indices one in front and one in the back (index1(R2) and index3(R3)). Goal is to de-multiplex the samples. 

Import modules like numpy and bioinfo
Define the reverse_complementary function() # the function doc itself is in answer_part2.md

At the same time open all the four files using Argparse or with open-read command:
    Loop over the lines in all the files one at a time and assign header, sequence, '+', and qscore to each line either using readline() command or the modulo function
    Once the file is open and read, you can either use the write out function or use the if/elif statement for organizing our reads into a)unknown&lowquality, b)Index hopping, c)Matched files)). More on this is coming up at the bottom
    Call the reverse complementary function to reverse and complement the DNA bases in the R3 or Index 2 file
    with open- write command write out all the  conditions stated
        Write out low quality file, matched files, and index hopped file: 
            if sequence in R2 and R3 has N and/or are mismatch at some base pairs (but not all) transfer that read into low quality file
            elif append the sequence in R2 and R3 to header line of R1 and R4 and write it in the matched file (2 files for Read1 and Read2)
                Open the two matched files at the same time and look for index hopping
                If Index1 and Index 2 are present in the 24 indexes and are not similar for a given read in R1 and R4. Pull/grab out that read from the matched file and now store it in the index hopped file. 
                Once you have removed index hopping now segregate matched files based on the 24 known indices.  
                By the end we will have 24+24+2(mismatch/lowquality)+2(index hopped files) = 52 files
                Also at the end write a loop to count the number of lines in matched index pair file, index hopped file and lowquality/mismatch file.