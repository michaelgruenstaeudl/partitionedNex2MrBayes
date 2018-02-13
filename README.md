*partitNex2MrBayes*: Multisequence alignments with charset definitions to analysis-ready MrBayes input
======================================================================================================

*partitNex2MrBayes* is a [Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) scrip to convert a multisequence DNA alignment in NEXUS format that contains character set definitions (hereafter 'charsets') to a partitioned MrBayes input file upon evaluating optimal partitioning schemes and best-fitting nucleotide substitution models.


Description of *partitNex2MrBayes*
----------------------------------

The script *partitNex2MrBayes* converts a multisequence DNA alignment in [NEXUS](https://en.wikipedia.org/wiki/Nexus_file) format that contains character set definitions (hereafter 'charsets') into a partitioned NEXUS file ready for analysis with [MrBayes](https://mrbayes.sourceforge.net). Alignment regions of different charsets are extracted, passed to one of several software tools for (a) evaluating optimal data partition schemes and (b) identifying best-fitting nucleotide substitution models. The alignment regions of different charsets are then concatenated again. Where necessary, the best-fitting nucleotide substitution models are converted to the specifications read by MrBayes. A command block for MrBayes is appended to the NEXUS file that integrates these partition definitions.

Software tools for evaluating optimal data partition schemes that *partitNex2MrBayes* can interact with:
- PartitionFinder2
- PartitionTest

Software tools for evaluating best-fitting nucleotide substitution models that *partitNex2MrBayes* can interact with:
- jModelTest
- SmartModelSelection
- Modeltest-NG


Details of the functioning of *partitNex2MrBayes*
----------------------------------

Several tests regarding the character set definitions are performed during the execution of the script. First, the script evaluates if the entire alignment is covered by partitions. If not, additional partitions are defined so that all partitions together form a continuous range from {1} to {total size of matrix}. Second, the script evaluates if any of the character definitions overlap with one another. If an overlap is detected, the script exists with an error. Consequently, the initial character set definitions of the NEXUS file do not need to cover the entire alignment length, but must not be overlapping.


Input/Output of *partitNex2MrBayes*
----------------------------------

## ARGS:
    input data file:   name of, and file path to, the input NEXUS file
    input config file: name of, and file path to, the input CONFIG file
    modeltesting type: name of modeltesting tool used (available: jmodeltest, partitionfinder, sms)
    modeltesting tool: name of, and file path to, the modeltesting binary or script
    output file:       name of, and file path to, the output file
## OUTP:
    MRBAYES file:      a NEXUS file ready for analysis with MrBayes



Status
------

Under active development
