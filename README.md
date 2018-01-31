# partitionedNex2MrBayes

Bash scrip to convert a partitioned NEXUS file to a partitioned MrBayes analysis file upon evaluating the best-fitting nucleotide substitution models

## Description
The script converts a DNA alignment in NEXUS format that contains character set ("charset") definitions into a partitioned NEXUS file ready for analysis with MrBayes. The character sets (hereafter "partitions") are hereby extracted, passed to one of several software tools for modeltesting (jModelTest, SmartModelSelection, PartitionFinder) to identify the best-fitting models of nucleotide substitition and then concatenated again. Where necessary, the best-fitting nucleotide substitution models are converted to the specifications read by MrBayes. A command block for MrBayes is appended to the NEXUS file that integrates these partition definitions.

Several tests regarding the character set definitions are performed during the execution of the script. First, the script evaluates if the entire alignment is covered by partitions. If not, additional partitions are defined so that all partitions together form a continuous range from {1} to {total size of matrix}. Second, the script evaluates if any of the character definitions overlap with one another. If an overlap is detected, the script exists with an error. Consequently, the initial character set definitions of the NEXUS file do not need to cover the entire alignment length, but must not be overlapping.

## Input and output
####ARGS:
    input data file:   name of, and file path to, the input NEXUS file
    input config file: name of, and file path to, the input CONFIG file
    modeltesting type: name of modeltesting tool used (available: jmodeltest, partitionfinder, sms)
    modeltesting tool: name of, and file path to, the modeltesting binary or script
    output file:       name of, and file path to, the output file
####OUTP:
    MRBAYES file:      a NEXUS file ready for analysis with MrBayes

## Status
In active development
