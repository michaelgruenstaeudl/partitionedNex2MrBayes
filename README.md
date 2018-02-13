*partitNex2MrBayes*
===================

A Bash script to convert a multisequence DNA alignment in NEXUS format with character set definitions to a partitioned MrBayes input file upon inferring optimal partitioning schemes and best-fitting nucleotide substitution models through third-party software tools.


Description of *partitNex2MrBayes*
----------------------------------

*partitNex2MrBayes* is a [Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) script to convert a  multisequence DNA alignment in [NEXUS](https://en.wikipedia.org/wiki/Nexus_file) format (hereafter 'matrix') that contains character set definitions (hereafter 'charsets') into a partitioned DNA matrix in NEXUS format ready for analysis with [MrBayes](https://mrbayes.sourceforge.net). The script hereby wraps and calls third-party software tools for inferring optimal data partitioning schemes and best-fitting nucleotide substitution models.

Software tools for inferring optimal data partitioning schemes wrapped by *partitNex2MrBayes*:
* PartitionFinder2 (link here)
* PartitionTest (link here)

Software tools for inferring best-fitting nucleotide substitution models wrapped by *partitNex2MrBayes*:
* jModelTest (link here)
* SmartModelSelection (link here)
* Modeltest-NG (link here)


Details of the analysis process of *partitNex2MrBayes*
------------------------------------------------------

### 1. Evaluation of presence and contiguity of charset definitions
Several evaluations regarding charset definitions are performed during the execution of the script. 
(a) The script evaluates if the entire alignment is covered by partitions. If not, additional partitions are defined so that all partitions together form a continuous range from {1} to {total size of matrix}.
(b) The script evaluates if any of the character definitions overlap with one another. If an overlap is detected, the script exists with an error. Consequently, the initial character set definitions of the NEXUS file do not need to cover the entire alignment length, but must not be overlapping.

### 2a. Using *partitNex2MrBayes* to infer optimal partitioning schemes
The script can be executed to call one of several software tools for inferring optimal partitioning schemes (i.e., PartitionFinder2, PartitionTest). Under this choice, the matrix of the input file is passed to the software tool selected in whole. Upon calling the software, the output is [complete here].

### 2b. Using *partitNex2MrBayes* to infer best-fitting nucleotide substitution models
The script can be executed to call one of several software tools for inferring best-fitting nucleotide substitution models (i.e., jModelTest, SmartModelSelection, Modeltest-NG, PartitionFinder2, PartitionTest). Under this choice, the matrix regions are extracted based on charset boundaries and then passed individually to the software tool selected. Where necessary (e.g., model selection via *PartitionTest* and *SmartModelSelection*), the best-fitting nucleotide substitution models inferred are converted to the precise model specifications read by MrBayes. Upon inference of best-fit models, matrix regions are concatenated again and saved as a reassembled NEXUS file.

### 3. Output is MrBayes-ready
Under each choice, a MrBayes command block is appended to the reassembled NEXUS file which specifies the partition and substitution model definitions. Thus, the output can be loaded into MrBayes without further adjustments.


Input/Output of *partitNex2MrBayes*
----------------------------------

### Input arguments
* __`input data file`__: name of, and file path to, the input NEXUS file
* __`input config file`__: name of, and file path to, the input CONFIG file
* __`modeltesting type`__: name of modeltesting tool used (available: jmodeltest, partitionfinder, sms)
* __`modeltesting tool`__: name of, and file path to, the modeltesting binary or script
* __`output file`__: name of, and file path to, the output file

### Output:
* __`MRBAYES file`__: a NEXUS file ready for analysis with MrBayes


Current Issues
--------------

Foo bar baz


Status
------

Under active development; pre-alpha
