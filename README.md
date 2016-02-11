# GND data analysis

This is to remind Jonathan what we're planning to do so he doesn't need to ask over and over!

23 animals, 1 control (#24), 4 repeats (half given treatment bifidobacteria)

4 extraction methods:
 - Using kit from post-enrichment
 - Boiled prep post-enrich
 - kit pre-enrichment
 - kit from faeces
 
 Controls: 4 synthetic libraries prepared from purified PCR products. Done using Boiled prep.
 
 Sequences associated with each serogroups - sequence
 
 gnd_seqs_metadata.xlsx (18thou counts is number of reads) librarynum_animalnum_S#counts
 
 1. Is there a relationship between e.coli diversity given treatment.
 2. Treatment effect?
 3. Sample process effect?
 
 2 thresholds of 15k, 18k. Repeat analysis with each and see if there's any difference.

 Can we cluster the serogroups into 'major' serogroups (essentially getting rid of the SLVs) and then look at the diversity stuff on that as well.
 
 TODOs:
 
 1. Read in data
    a. Merge the .fa/.txt sheets up to give an id/counts data frame, and id/sequence data frame.
    b. Merge the metadata.xlsx sheet into the id/counts data frame
 2. Compute pair-wise distances between serogroups (the sequences).
 3. Use that to look at clustering the serogroups and see what we get.
    a. Is there a small number of groups with lots of SLVs?
    b. Do we get differences if we include the counts per sequence across all the cows (to get some sort of 'importance' measure?)

 4. Compute pair-wise distances between samples (23*4)? Not sure this is exactly what we want, or whether we need only the differences between the reads (thus sequences).
 5. Use PERMANOVA to do testing there. If we're on sequences, that could be fun, as we have a few thousand reads per sample, so 433000 rows. May need our own PERMANOVA routine here,
 as the distance matrix is sparse, in that it is blocked (there's only 1070x1070 bits of info). Also in the PERMANOVA we're only interesting in shuffling samples in order to test
 whether method makes sense or not.
 6. Repeat with both 15k and 18k. Is there a difference?
