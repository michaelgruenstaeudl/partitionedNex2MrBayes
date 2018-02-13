*partitNex2MrBayes*
===================

A Bash script to convert a multisequence DNA alignment in NEXUS format with character set definitions to a partitioned MrBayes input file upon inferring optimal partitioning schemes and best-fitting nucleotide substitution models through third-party software tools.


Description of *partitNex2MrBayes*
----------------------------------

*partitNex2MrBayes* is a [Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) script to convert a  multisequence DNA alignment in [NEXUS](https://en.wikipedia.org/wiki/Nexus_file) format (hereafter 'matrix') that contains character set definitions (hereafter 'charsets') into a partitioned DNA matrix in NEXUS format ready for analysis with [MrBayes](https://mrbayes.sourceforge.net). The script hereby wraps and calls third-party software tools for inferring optimal data partitioning schemes and best-fitting nucleotide substitution models.

Software tools for inferring optimal data partitioning schemes wrapped by *partitNex2MrBayes*:
* PartitionFinder2 (https://github.com/brettc/partitionfinder)
* PartitionTest (https://github.com/ddarriba/partitiontest)

Software tools for inferring best-fitting nucleotide substitution models wrapped by *partitNex2MrBayes*:
* jModelTest (https://github.com/ddarriba/jmodeltest2)
* SmartModelSelection (http://www.atgc-montpellier.fr/sms/binaries.php)
* Modeltest-NG (https://github.com/ddarriba/modeltest)


Analysis process of *partitNex2MrBayes*
---------------------------------------

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

Control of *partitNex2MrBayes*
------------------------------

Analyses conducted via *partitNex2MrBayes* are controlled by
(a) commandline parameters during script execution, and
(b) a configuration file (default: *partitNex2MrBayes.cfg*).

The configuration file controls the settings of the third-party software tools that *partitNex2MrBayes* wraps and calls for analysis. These settings include the choice of optimality criterion (e.g., AIC, AICc, BIC) and the choice of search algorithm (e.g., greedy, heuristic clustering) and are largely idiosyncratic to the wrapped software tool.


Input/Output of *partitNex2MrBayes*
----------------------------------

### Input

#### Mandatory commandline arguments
* __`input data file`__ (commandline option __`-f`__): name of, and file path to, the input NEXUS file. No default exists.
* __`input config file`__ (commandline option __`-c`__): name of, and file path to, the input CONFIG file. No default exists.
* __`modeltesting type`__ (commandline option __`-t`__): name of partitionfinding/modeltesting tool selected. Available: jmodeltest, modeltest_ng, partitionfinder, partitiontest, sms. No default exists.
* __`modeltesting tool`__ (commandline option __`-b`__): name of, and file path to, the partitionfinding/modeltesting binary or script. No default exists.
* __`output file`__ (commandline option __`-o`__): name of, and file path to, the output file. No default exists.

#### Optional commandline switches
* commandline switch __`-u`__: Estimating model fit only for the user-defined partitioning scheme. Only applicable for PartitionFinder. Default is off.
* commandline switch __`-k`__: Keeping temporary files. Default is off.
* commandline switch __`-v`__: Displays verbose status messages. Default is off.


### Output
* __`MRBAYES file`__: an output NEXUS file ready for analysis with MrBayes. The file comprises a DATA block, a commented-out SETS Block and a MRBAYES block.


Example Usage
-------------

##### 1. Setting input file and working directory
``` 
INFILE=/home/user/analysis/Gruenstaeudl_2017.nex
INFSTEM=${INFILE%.nex*}
WORKDIR=/home/user/analysis
```

##### 2. Setting script directory and script file
```
PN2MB_SH=/home/user/git/partitNex2MrBayes/partitNex2MrBayes.sh
PN2MB_CFG=/home/user/git/partitNex2MrBayes/partitNex2MrBayes.cfg
```

##### 3. Setting third-party software tool
```
MDLTSTTYPE=partitionfinder
PATH_TO_PARTITIONFINDER=/home/user/git/partitionfinder/PartitionFinder.py
OUTSTEM=${INFSTEM}_${MDLTSTTYPE}
```

##### 4. Execute *partitNex2MrBayes*
```
bash $PN2MB_SH -f $INFILE -c $PN2MB_CFG -t $MDLTSTTYPE -b $PATH_TO_PARTITIONFINDER -o $WORKDIR/${OUTSTEM}.mrbayes -v -k > $WORKDIR/${OUTSTEM}.log
```

Current Issues
--------------
* 2017-11: Modeltest-NG crashes at times. The script *partitNex2MrBayes* does not recognize such a crash and waits indefinitely. (Note to self: Solve this issue by having a maximum execution time and by selecting the general GTR+I+G model by default in case of such a crash.) 


Status
------

Under active development; pre-alpha
