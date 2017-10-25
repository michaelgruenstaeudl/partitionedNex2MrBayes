#!/bin/bash

########################################################################

THIS_SCRIPT=`basename ${BASH_SOURCE[0]}`
#DESCRIPTION="A bash shell scrip to convert a partitioned NEXUS file into a partitioned MrBayes analysis file"
AUTHOR="Michael Gruenstaeudl, PhD"
#COPYRIGHT="Copyright (C) 2015-2017 $AUTHOR"
#CONTACT="m.gruenstaeudl@fu-berlin.de"
VERSION="2017.10.23.1400"
USAGE="bash $THIS_SCRIPT -f <path_to_NEXUS_file> -c <path_to_config_file> -t <type_of_modeltesting_tool> -b <path_to_modeltesting_tool> -o <path_to_outfile>"

########################################################################

#    The script converts a DNA alignment in NEXUS format that contains character set ("charset") definitions into a partitioned NEXUS file ready for analysis with MrBayes. The character sets (hereafter "partitions") are hereby extracted, passed to one of several software tools for modeltesting (jModelTest, SmartModelSelection, PartitionFinder) to identify the best-fitting models of nucleotide substitition and then concatenated again. Where necessary, the best-fitting nucleotide substitution models are converted to the specifications read by MrBayes. A command block for MrBayes is appended to the NEXUS file that integrates these partition definitions.
#    Several tests regarding the character set definitions are performed during the execution of the script. First, the script evaluates if the entire alignment is covered by partitions. If not, additional partitions are defined so that all partitions together form a continuous range from {1} to {total size of matrix}. Second, the script evaluates if any of the character definitions overlap with one another. If an overlap is detected, the script exists with an error. Consequently, the initial character set definitions of the NEXUS file do not need to cover the entire alignment length, but must not be overlapping.
#    ARGS:
#        input data file:   name of, and file path to, the input NEXUS file
#        input config file: name of, and file path to, the input CONFIG file
#        modeltesting type: name of modeltesting tool used (available: jmodeltest, partitionfinder, sms)
#        modeltesting tool: name of, and file path to, the modeltesting binary or script
#        output file:       name of, and file path to, the output file
#    OUTP:
#        MRBAYES file:      a NEXUS file ready for analysis with MrBayes

########################################################################

# TO DO:
# [currently nothing]

# PAPER: BMC Evol Biol. 2010 Aug 9;10:242. doi: 10.1186/1471-2148-10-242. Performance of criteria for selecting evolutionary models in phylogenetics: a comprehensive study based on simulated datasets. Luo A1, Qiao H, Zhang Y, Shi W, Ho SY, Xu W, Zhang A, Zhu C.
# TELLS THAT BIC SHOULD BE USED FOR MODEL SELECTION.
# 

########################################################################

# INITIAL SETTINGS AND COMMANDLINE INTERFACE

# Toggle on: throw errors
set -e

# Start timing execution time
RUNTIME_start=$(date +%s)
RUNTIME_start_pp=$(date '+%Y.%m.%d.%H%M')

# Initialize variables to default values
OPT_A="name of, and file path to, the input NEXUS file"
OPT_B="name of, and file path to, the input config file"
OPT_C="name of modeltesting software tool (available: jmodeltest, partitionfinder, sms)"
OPT_D="name of, and file path to, the modeltesting binary or script"
OPT_E="name of, and file path to, the output file"
keepflsBOOL=0
verboseBOOL=0
usrschmBOOL=0

# Set fonts for help function
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# Help function
function my_help {
    echo -e \\n"Help documentation for ${BOLD}$THIS_SCRIPT${NORM}."\\n
    echo -e "Version: $VERSION | Author: $AUTHOR"\\n
    echo -e "${REV}Usage:${NORM} $USAGE"\\n
    echo "MANDATORY command line switches:"
    echo "${REV}-f${NORM}  --Sets name of, and file path to, ${BOLD}NEXUS file${NORM}. No default exists."
    echo "${REV}-c${NORM}  --Sets name of, and file path to, ${BOLD}config file${NORM}. No default exists."
    echo "${REV}-t${NORM}  --Sets ${BOLD}type${NORM} of modeltesting software tool. Available: jmodeltest, partitionfinder, sms. No default exists."
    echo "${REV}-b${NORM}  --Sets name of, and file path to, the modeltesting ${BOLD}binary or script${NORM}. No default exists."
    echo "${REV}-o${NORM}  --Sets name of, and file path to, ${BOLD}output file${NORM}. No default exists."
    echo ""
    echo "OPTIONAL command line switches:"
    echo "${REV}-u${NORM}  --Estimating model fit only for the ${BOLD}user-defined${NORM} partitioning scheme. Only applicable for PartitionFinder. Default is ${BOLD}off${NORM}."
    echo "${REV}-k${NORM}  --${BOLD}Keeps${NORM} temporary files. Default is ${BOLD}off${NORM}."
    echo "${REV}-v${NORM}  --Displays ${BOLD}verbose${NORM} status messages. Default is ${BOLD}off${NORM}."
    echo -e "${REV}-h${NORM}  --Displays this help message."\\n
    exit 1
}

# Check the number of arguments. If none are passed, print help and exit.
numArgmts=$#
if [ $numArgmts -eq 0 ]; then
    my_help
fi

########################################################################

# GETOPTS CODE

while getopts :f:c:t:b:o:hukv FLAG; do
  case $FLAG in
    f)  OPT_A=$OPTARG
        ;;
    c)  OPT_B=$OPTARG
        ;;
    t)  OPT_C=$OPTARG
        ;;
    b)  OPT_D=$OPTARG
        ;;
    o)  OPT_E=$OPTARG
        ;;
    h)  my_help
        ;;
    u)  usrschmBOOL=1
        ;;
    k)  keepflsBOOL=1
        ;;
    v)  verboseBOOL=1
        ;;
    \?) echo 'Invalid option: -$OPTARG' >&2
        exit 1
        ;;
    :)  echo 'Option -$OPTARG requires an argument.' >&2
        exit 1
        ;;
  esac
done

#shift $((OPTIND-1))

########################################################################
########################################################################

# GENERAL FUNCTIONS

get_current_time() {
    date '+%Y-%m-%d %Z %H:%M:%S'
}

assemble_mrbayes_block_general()
#   This function appends a MrBayes block to a file ($NamOTFILE).
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of outfile ($NamOTFILE)
#         $3: name of file with lset definitions ($lsetDefns)
#         $4: name of file with charsets lines ($charSetsL)
#         $5: name of config file ($reformCFG)
#   OUP:  none; operates on lset definitions file
{
    # INTERNAL FUNCTION CHECKS
    (($# == 5)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    charsetLines=$(cat $1/$4)
    echo -e 'begin mrbayes;' >> $2
    echo -e 'set autoclose=yes;' >> $2
    echo -e 'set nowarnings=yes;' >> $2
    echo "$charsetLines" >> $2
    echo -n 'partition combined =' $(echo "$charsetLines" | wc -l) ': ' >> $2 # Line must not end with linebreak, thus -n
    echo "$charsetLines" | awk '{print $2}' | awk '{ORS=", "; print; }' >> $2
    sed -i 's/\(.*\)\,/\1;/' $2 # Replace final comma with semi-colon; don't make this replacement global
    echo -e '\nset partition = combined;' >> $2 # Option -e means that \n generates new line
    cat $1/$3 >> $2
    echo 'prset applyto=(all) ratepr=variable;' >> $2
    echo 'unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);' >> $2
    echo -e 'end;\n' >> $2
}

check_inp_outp_availability()
#   This function checks the availability of input and output files.
#   INP:  $1: path and name of input: NEXUS file ($nexINFILE)
#         $2: path and name of input: CONFIG file ($cfgINFILE)
#         $3: path and name of input: modeltest binary/script ($BinMDLTST)
#         $4: path and name of output file ($NamOTFILE)
#   OUP:  none
{
    # INTERNAL FUNCTION CHECKS
    (($# == 4)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # Checking if input files exists
    #if [[ ! -f $1 ]]; then 
    if [ ! -f $1 ]; then 
        echo -ne " ERROR | $(get_current_time) | NEXUS file not found: $1\n"
        exit 1
    fi
    if [ ! -f $2 ]; then 
        echo -ne " ERROR | $(get_current_time) | Config file not found: $2\n"
        exit 1
    fi
    if [ ! -f $3 ]; then 
        echo -ne " ERROR | $(get_current_time) | Modeltest binary/script file not found: $3\n"
        exit 1
    fi
    # Check if outfile already exists
    #if [[ $3 = /* ]]; then # File is generated by GETOPTS, so checking its presence is not useful
    if [ -f $4 ]; then
        echo -ne " ERROR | $(get_current_time) | Outfile already exists in filepath: $4\n"
        exit 1
    fi
}

ensure_partitions_form_continuous_range()
#   This function tests if the partitions form a continuous range. 
#   If the partitions don't form such a continuous range, additional
#   ranges are inserted such that all ranges taken together form a 
#   continuous range from {1} to {total size of matrix}.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file containing the SETS-block
#         $3: total size of matrix
#   OUP:  update of $1
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # Get charset definition lines
    charsetLines=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $1/$2 | grep 'charset')
    # Convert from discrete to continuous range
    continuousRange1=$(echo "$charsetLines" | awk '{ split($4,curr,/[-;]/); currStart=curr[1]; currEnd=curr[2] } currStart > (prevEnd+1) { print "charset " "new"++cnt " = " prevEnd+1 "-" currStart-1 ";" } { print; prevEnd=currEnd }')
    # Add concluding partition, if missing
    matrxLngth=$3
    continuousRange2=$(echo "$continuousRange1" | tail -n1 | awk -F'[[:space:]]*|-' -v lngth=$matrxLngth ' $5<lngth {print "charset newFinal = " $5+1 "-" lngth ";"}')
    # Update SETS-block
    echo "begin sets;" > $1/$2
    echo "$continuousRange1" >> $1/$2 # NOTE: $continuousRange1 is a multi-line variable, whose linebreaks shall be maintained; thus, it is passed as doublequote-enclosed.
    echo "$continuousRange2" >> $1/$2
    #echo "end;" >> $1/$2
    echo -n "end;" >> $1/$2
}

get_num_of_avail_cores()
#   This function evaluates the numberof available cores on a system
#   INP:  none
#   OUP:  returns integer
{
    nmbrCores=$(nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null)  # NOTE: The "||" indicates that if "nproc" fails, do "grep ^processor /proc/cpuinfo".
    posIntgr='^[0-9]+$'
    if ! [[ $nmbrCores =~ $posIntgr ]]; then
        nmbrCores=1
    fi
    echo $nmbrCores
}

reconstitute_datablock_as_interleaved()
#   This function reconstitutes a data block from individual partitions, but makes the data block interleaved.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of DATA-block ($dataBlock)
#         $3: name of outfile ($NamOTFILE)
#   OUP:  none; writes to outfile
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # Step 1: Append dimension and format info of DATA-block
    cat $1/$2 | sed -n '/begin data\;/,/matrix/p' >> $3
    # Step 2: Specify that matrix is interleaved
    sed -i '/format/ s/\;/ interleave\=yes\;/' $3 # in-line replacement of semi-colon with interleave-info plus semi-colon, if line has keyword 'format'
    # Step 3: Append the individual partitions
    for p in $(ls $1/partition_[0-9]* | grep -v "\.bestModel\|\.phy\|\.runlog\|_sms\.csv"); do # NOTE: Use [0-9] to separate from 'partition_finder.cfg'
        echo -ne "\n[$(basename $p)]\n" >> $3
        pureMatrx=$(cat $p | sed -n '/matrix/{:a;n;/;/b;p;ba}')
        algnMatrx=$(echo "$pureMatrx" | column -t)
        echo "$algnMatrx" >> $3 # Append only the matrix of a partition, not the preceeding dimension and format info; also, don't append the closing ';\nend;' per partition
    done
    # Step 4: Append a closing ';\nend;'
    echo -e "\n;\nend;" >> $3
}

reformat_config_file()
#   This function performs the following steps:
#   (a) it removes comments from a CONFIG file,
#   (b) it removes blank lines (both blank by whitespaces and by tabs),
#   (c) it removes the configurations of modeltesting types not selected.
#   INP:  $1: name of complete CONFIG file ($cfgINFILE)
#         $2: type of modeltesting selected: ($typMDLTST)
#         $3: name of a reformatted CONFIG file ($tmpFOLDER/$reformCFG)
#   OUP:  the re-formatted CONFIG file
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    ## REMOVES LEADING WHITESPACES FROM EVERY LINE
    #  (so that formatting (i.e., indets) does not matter)
    configWithoutLeadWhitesp=$(sed -e 's/^[ \t]*//' $1)
    
    ## REMOVES COMMENTS FROM A CONFIG FILE
    #  (i.e., delete all info that starts with a pound sign)
    configWithoutComments=$(echo "$configWithoutLeadWhitesp" | grep -v "^#")

    ## REMOVES BLANK LINES
    #  (both blank by whitespaces and by tabs)
    configWithoutBlankLines=$(echo "$configWithoutComments" | sed '/^\s*$/d')
    
    ## KEEP ONLY CONFIGURATION OF MODELTEST TYPE SELECTED
    configOnlyTypeSelected=$(echo "$configWithoutBlankLines" | grep -A1 --ignore-case --no-group-separator "\%$2\|\%mrbayes")
    
    echo "$configOnlyTypeSelected" > $3
}


reformat_nexus_file()
#   This function performs the following steps:
#   (a) it removes comments from a NEXUS file,
#   (b) it removes blank lines (both blank by whitespaces and by tabs), 
#   (c) it standardizes critical command lines (i.e. begin data, 
#       begin sets, dimensions, etc.) as lowercase, and
#   (d) it removes lines from the sets-block that do not start with the keyword "charset"
#   INP:  $1: name of complete NEXUS file ($nexINFILE)
#         $2: name of a reformatted NEXUS file ($tmpFOLDER/$reformNEX)
#   OUP:  the re-formatted NEXUS file
{
    # INTERNAL FUNCTION CHECKS
    (($# == 2)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    ## REMOVES COMMENTS FROM A NEXUS FILE
    #  (i.e., delete all info enclosed by square brackets)
    nexusWithoutComments=$(sed 's/\[[^]]*\]//g' $1)

    ## REMOVES BLANK LINES
    #  (both blank by whitespaces and by tabs)
    nexusWithoutBlankLines=$(echo "$nexusWithoutComments" | sed '/^\s*$/d')

    ## REMOVES LEADING WHITESPACES FROM EVERY LINE
    #  (so that formatting (i.e., indets) does not matter)
    nexusWithoutLeadWhitesp=$(echo "$nexusWithoutBlankLines" | sed -e 's/^[ \t]*//')

    ## STANDARDIZES CRITICAL COMMAND LINES AS LOWERCASE (i.e., converts 
    ## ENTIRE LINE that starts with a keyword TO LOWERCASE)
    #  NOTE: These conversions to lowercase are critical for the correct 
    #        section identifikation below.
    #  keyword1: "begin", line ends with semi-colon
    #  NOTE: the following command converts "Begin Data;" to "begin data;", not just "begin"!
    reformatForKeyword1=$(echo "$nexusWithoutLeadWhitesp" | awk 'BEGIN{IGNORECASE=1} /^ *begin\>.*; *$/ {$0=tolower($0)} 1')
    #  keyword2: "end;"
    reformatForKeyword2=$(echo "$reformatForKeyword1" | awk 'BEGIN{IGNORECASE=1} /^ *end\;/ {$0=tolower($0)} 1')
    #  keyword3: "matrix"
    reformatForKeyword3=$(echo "$reformatForKeyword2" | awk 'BEGIN{IGNORECASE=1} /^ *matrix\>/ {$0=tolower($0)} 1')
    #  keyword4: "dimensions", line ends with semi-colon
    reformatForKeyword4=$(echo "$reformatForKeyword3" | awk 'BEGIN{IGNORECASE=1} /^ *dimensions\>.*; *$/ {$0=tolower($0)} 1')
    #  keyword5: "format", line ends with semi-colon
    reformatForKeyword5=$(echo "$reformatForKeyword4" | awk 'BEGIN{IGNORECASE=1} /^ *format\>.*; *$/ {$0=tolower($0)} 1')
    #  Convert only specific keyword to lowercase
    reformatForKeyword6=$(echo "$reformatForKeyword5" | awk 'tolower($1)=="charset"{$1=tolower($1)}1')

    ## SORT CHARSETS
    charsetLines=$(echo "$reformatForKeyword6" | sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' | grep 'charset')
    charsetLinesSorted=$(echo "$charsetLines" | sort -n -t'=')
    
    ## REMOVES LINES FROM THE SETS-BLOCK THAT DO NOT START WITH THE KEYWORD "CHARSET"
    #  NOTE: This check must come after the conversion to lowercase chars!
    #  NOTE: This step generates the outfile (i.e., reformatted NEXUS file)
    echo -e "#NEXUS\n" > $2
    echo "$reformatForKeyword6" | sed -n '/begin data\;/,/end\;/p' >> $2
    echo -e "\nbegin sets;" >> $2
    echo "$charsetLinesSorted" >> $2
    echo "end;" >> $2

    ## CHECK THAT DATATYPE IS "DNA|dna"
    #  NOTE: This check must come after the conversion to lowercase chars!
    dataTypeB=$(cat $2 | grep format | sed -n 's/.*datatype=\([^ ]*\).*/\1/p') # Extract info on datatype
    if [ ! "$dataTypeB" = "dna" ]; then
        echo -ne " ERROR | $(get_current_time) | Datatype not set as: DNA\n"
        exit 1
    fi
}

split_matrix_into_partitions()
#   This function takes the matrix of a NEXUS file and extracts from it a 
#   character range (hereafter "sub-matrix") specified by the definition 
#   of a charsets line.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file containing the DATA-block
#         $3: name of file containing the SETS-block
#   OUP:  file with name $charrngFn in temp folder
#         file with name $partFname in temp folder
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    pureMatrx=$(sed -n '/matrix/{:a;n;/;/b;p;ba}' $1/$2)
    charsetLines=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $1/$3 | grep 'charset')
    nexusNew1='#NEXUS\n\n'
    nexusNew2=$(sed -e '/matrix/,$d' $1/$2) # Get section between "#NEXUS" and "MATRIX"
    nexusNew3='\nmatrix'
    nexusNew4=';\nend;\n'
    myCounter=0
    while IFS= read -r line; do 
        myCounter=$((myCounter+1))
        partitFn1=$(printf %03d $myCounter)
        partitFn2=$(echo "$line" | awk '{print $2}')
        partitnFn=partition_${partitFn1}_${partitFn2}
        # Get the info on the charset range
        charsetRanges=$(echo "$line" | awk '{print $4}' | sed 's/\;//')
        # Step 1 of assembling the new partition file
        echo -e "$nexusNew1$nexusNew2$nexusNew3" > $1/$partitnFn
        # Step 2 of assembling the new partition file: Add the sub-matrix
        awk 'NR==FNR{start=$1;lgth=$2-$1+1;next} {print $1, substr($2,start,lgth)}' FS='-' <(echo "$charsetRanges") FS=' ' <(echo "$pureMatrx") >> $1/$partitnFn
        # Step 3 of assembling the new partition file
        echo -e "$nexusNew4" >> $1/$partitnFn
        # Get the length of the sub-matrix
        mtrxLngth=$(awk '{print $2-$1+1}' FS='-' <(echo "$charsetRanges"))
        # Replace the number of characters with the length of the sub-matrix
        sed -i "/dimensions / s/nchar\=.*\;/nchar\=$mtrxLngth\;/" $1/$partitnFn
    done <<< "$charsetLines" # Using a here-string
}

split_nexus_into_blocks()
#   This function splits a NEXUS file into individual blocks.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of reformatted NEXUS file ($reformNEX)
#         $3: name of file containing the extracted DATA-block ($dataBlock)
#         $4: name of file containing the extracted SETS-block ($setsBlock)
#   OUP:  file with DATA-block
#         file with SETS-block
{
    # INTERNAL FUNCTION CHECKS
    (($# == 4)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # EXTRACT THE BLOCKS
    cat $1/$2 | sed -n '/begin data\;/,/end\;/p' > $1/$3
    cat $1/$2 | sed -n '/begin sets\;/,/end\;/p' > $1/$4
    # CHECK IF DATABLOCK SUCCESSFULLY GENERATED
    if [ ! -s "$1/$3" ]; then
        echo -ne " ERROR | $(get_current_time) | No file size: $1/$3\n"
        exit 1
    fi
    # CHECK IF SETSBLOCK SUCCESSFULLY GENERATED
    if [ ! -s "$1/$4" ]; then
        echo -ne " ERROR | $(get_current_time) | No file size: $1/$4\n"
        exit 1
    fi
}

test_if_partitions_overlap()
#   This function tests if any of the partitions (i.e., charset ranges)
#   overlap. If they do, an error is thrown.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file containing the SETS-block
#   OUP:  none
{
    # INTERNAL FUNCTION CHECKS
    (($# == 2)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    charsetLines=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $1/$2 | grep 'charset')
    charsetRanges=$(echo "$charsetLines" | awk '{print $4}' | sed 's/\;//')
    # Test if any of the charset ranges are overlapping
    charsetOverlap=$(echo "$charsetRanges" | awk -F"-" 'Q>=$1 && Q{print val}{Q=$NF;val=$0}')
    if [ ! "$charsetOverlap" = "" ]; then 
        echo -ne " ERROR | $(get_current_time) | Charset range overlaps with subsequent range: $charsetOverlap\n"
        exit 1
    fi
}

########################################################################
########################################################################

# FUNCTIONS SPECIFIC TO 'JMODELTEST'

convert_models_into_lset()
#   This function converts the names of nucleotide substitution models 
#   into lset definitions readable by MrBayes.
#   NOTE: It also adds the model names as comments to the end of each 
#   lset line to allow the applyto-assignment (at the end of the function)!
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file with best models ($bestModls)
#         $3: name of file with lset definitions file ($lsetDefns)
#         $4: name of file with charsets lines ($charSetsL)
#   OUP:  none; operates on lset definitions file
#   NOTE:   The model abbreviation (e.g., `GTR+I+G`) must be at the line 
#           end and no whitespace must come thereafter. This end of line 
#           position precludes that string `GTR+I+G` would be replaced by 
#           `GTR+I`.
{
    # INTERNAL FUNCTION CHECKS
    (($# == 4)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    if [ ! -f $1/$2 ]; then 
        echo -ne " ERROR | $(get_current_time) | File not found: $1/$2\n"
        exit 1
    fi

    # STEP 1: Convert the model names to lset definitions
    while read LINE
    do
        modelKeyw=$(echo "$LINE" | awk '{print $2}' | cut -c1-2)
        if [ "$modelKeyw" == "GT" ]; then #GTR
            # NOTE: The dollar-sign behind the model name is critical to avoid that, for example, GTR+I replaces GTR+I+G
            # NOTE: The final ampersand (in \[\#&\]) copies the model name to the end of the line; this is important to allow the applyto-assignment!
            LINE=$(echo "$LINE" | sed 's/.* GTR$/lset applyto=() nst\=6\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* GTR+I$/lset applyto=() nst\=6 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* GTR+G$/lset applyto=() nst\=6 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* GTR+I+G$/lset applyto=() nst\=6 rates\=invgamma\; \[\#&\]/')
        elif [ "$modelKeyw" == "HK" ]; then #HKY
            LINE=$(echo "$LINE" | sed 's/.* HKY$/lset applyto=() nst\=2\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY+I$/lset applyto=() nst\=2 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY+G$/lset applyto=() nst\=2 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY+I+G$/lset applyto=() nst\=2 rates\=invgamma\; \[\#&\]/')
            # In SMS, "HKY" is called "HKY85"
            LINE=$(echo "$LINE" | sed 's/.* HKY85$/lset applyto=() nst\=2\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY85+I$/lset applyto=() nst\=2 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY85+G$/lset applyto=() nst\=2 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | sed 's/.* HKY85+I+G$/lset applyto=() nst\=2 rates\=invgamma\; \[\#&\]/')
        elif [ "$modelKeyw" == "K8" ]; then #K80(=K2P)
            LINE=$(echo "$LINE" | awk '/K80$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* K80$/lset applyto=() nst\=2\; \[\#&\]/') # NOTE: `/prset/!` means conduct replacement unless keyword ´prset´ present in line
            LINE=$(echo "$LINE" | awk '/K80\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* K80+I$/lset applyto=() nst\=2 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/K80\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* K80+G$/lset applyto=() nst\=2 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/K80\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* K80+I+G$/lset applyto=() nst\=2 rates\=invgamma\; \[\#&\]/')
        elif [ "$modelKeyw" == "SY" ]; then #SYM
            LINE=$(echo "$LINE" | awk '/SYM$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* SYM$/lset applyto=() nst\=6\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/SYM\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* SYM+I$/lset applyto=() nst\=6 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/SYM\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* SYM+G$/lset applyto=() nst\=6 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/SYM\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' "$LINE")
            LINE=$(echo "$LINE" | sed '/prset/! s/.* SYM+I+G$/lset applyto=() nst\=6 rates\=invgamma\; \[\#&\]/')
        elif [ "$modelKeyw" == "F8" ]; then #F81
            LINE=$(echo "$LINE" | sed 's/.* F81$/lset applyto=() nst\=1\; \[\#&\]/' "$LINE")
            LINE=$(echo "$LINE" | sed 's/.* F81+I$/lset applyto=() nst\=1 rates\=propinv\; \[\#&\]/' "$LINE")
            LINE=$(echo "$LINE" | sed 's/.* F81+G$/lset applyto=() nst\=1 rates\=gamma\; \[\#&\]/' "$LINE")
            LINE=$(echo "$LINE" | sed 's/.* F81+I+G$/lset applyto=() nst\=1 rates\=invgamma\; \[\#&\]/' "$LINE")
        elif [ "$modelKeyw" == "JC" ]; then #JC
            LINE=$(echo "$LINE" | awk '/JC\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* JC+I+G$/lset applyto=() nst\=1 rates\=invgamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/JC\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* JC+G$/lset applyto=() nst\=1 rates\=gamma\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/JC\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* JC+I$/lset applyto=() nst\=1 rates\=propinv\; \[\#&\]/')
            LINE=$(echo "$LINE" | awk '/JC$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1')
            LINE=$(echo "$LINE" | sed '/prset/! s/.* JC$/lset applyto=() nst\=1\; \[\#&\]/')
        else # In case a modelname is not recognized (e.g., TN93 with SMS), substitute with the default GTR+I+G
            LINE=$(echo "$LINE" | sed 's/.* GTR+I+G$/lset applyto=() nst\=6 rates\=invgamma\; \[\#& \- model not recognized \- replaced with GTR\+I\+G \]/')
        fi
        echo "$LINE" >> $1/$3
    done < $1/$2

    # STEP 2: Replace the partition names in the lset definitions with the line numbers in the charsets
    # NOTE: THIS IS WHY ORDERED PARTITIONS ARE CRITICAL!
    charsetLines=$(cat $1/$4)
    for charset in $(echo "$charsetLines" | awk '{print $2}'); do # Use doublequote-enclosing around a variable to maintain its newlines!
        lineNumbr=$(echo "$charsetLines" | awk '/'$charset'/{print NR; exit}')
        sed -i "/$charset/ s/applyto\=()/applyto\=($lineNumbr)/" $1/$3 # Double quotes are critical here due to variable
    done
}

extract_modelinfo_from_jmodeltest()
#   This function extracts modelinfo from the output of jModeltest; the output was saved into logfiles ending with the name ".bestModel"
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file to collect bestfit model per partition ($bestModls)
#         $3: boolean variable if verbose or not ($verboseBOOL)
#   OUP:  writes output to model overview file
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    count=$(ls -1 $1/*.bestModel 2>/dev/null | wc -l)
    if [ $count != 0 ]; then
        for file in ./$1/*.bestModel; do 
            echo -ne "$(basename $file)" | sed 's/partition_//g' | sed 's/\.bestModel//g' >> $1/$2
            bestModel=$(cat "$file" | grep -A1 ' Model selected:' | tail -n1 | sed -n 's/.*Model \= \([^ ]*\).*/\1/p' 2>/dev/null) # Parse best model from jModeltest output
            if [ ! -z "$bestModel" ]; then
                echo " $bestModel" >> $1/$2  # NOTE: model has to be preceeded by one space!
            else
                echo " GTR+I+G" >> $1/$2  # If error in modeltesting or parsing, set GTR+I+G as default model
                if [ $3 -eq 1 ]; then
                    echo -ne " WARN  | $(get_current_time) | No model extracted from file $file; setting model GTR+I+G as replacement.\n"
                fi
            fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No result files of the modeltesting (*.bestModel) found.\n"
        exit 1
    fi
}

modeltesting_via_jmodeltest()
#   This function conducts modeltesting via jModeltest for a series of input files; these input files start with the name "partition_"
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of jModeltest binary ($BinMDLTST)
#         $3: name of config file ($reformCFG)
#         $4: boolean variable if verbose or not ($verboseBOOL)
#   OUP:  none; generates *.bestModel* output files in temp folder
{
    # INTERNAL FUNCTION CHECKS
    (($# == 4)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done

    # USING CONFIG TEMPLATE
    CMDline=$(grep -A1 --ignore-case '^%JMODELTEST:' $1/$3 | tail -n1)
    CMDline=${CMDline//\[PATH_TO_JMODELTEST_JAR\]/$2} # Substitution with Bash's built-in parameter expansion
    CMDline=$(eval echo $CMDline)
    # MODELTESTING
    count=$(ls -1 ./$1/partition_* 2>/dev/null | wc -l)
    if [ $count != 0 ]; then
        for partit in ./$1/partition_*; do
        partit_RUNTIME_start=$(date +%s)
        coreCMD=${CMDline//\[PARTITION_NAME\]/$partit} # Substitution with Bash's built-in parameter expansion; # NOTE: Don't overwrite the variable $CMDline to perform different substitutions in each loop, instead save it as a new variable ($coreCMD)
        eval "$coreCMD 1>${partit}.bestModel 2>&1" # Executin of main_cmd
        partit_RUNTIME_dur=$(bc -l <<< "($(date +%s)-$partit_RUNTIME_start)/60")
        LC_ALL=C partit_RUNTIME_dur_pp=$(printf "%.3f minutes\n" $partit_RUNTIME_dur)
        if [ $4 -eq 1 ]; then
            echo -ne " INFO  | $(get_current_time) |   Processing complete for partition: $(basename $partit); Execution time: $partit_RUNTIME_dur_pp\n"
        fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No partition files (partition_*) found.\n"
        exit 1
    fi
}

########################################################################
########################################################################

# FUNCTIONS SPECIFIC TO 'SMS'

extract_modelinfo_from_sms()
#   This function extracts modelinfo from the output of sms; the output was saved into a file called "sms.csv" (default behaviour of SMS)
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file to collect bestfit model per partition ($bestModls)
#         $3: boolean variable if verbose or not ($verboseBOOL)
#   OUP:  writes output to model overview file
{
    # INTERNAL FUNCTION CHECKS
    (($# == 3)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    count=$(ls -1 $1/*_sms.csv 2>/dev/null | wc -l)
    if [ $count != 0 ]; then
        for file in ./$1/*_sms.csv; do 
            echo -ne "$(basename $file)" | sed 's/partition_//g' | sed 's/\_sms\.csv//g' >> $1/$2
            bestModel=$(cat "$file" | grep -A1 'Model;Decoration' | tail -n1 | awk -F';' '{print $1$2}' 2>/dev/null) # Parse best model from SMS output
            if [ ! -z "$bestModel" ]; then
                echo " $bestModel" >> $1/$2  # NOTE: model has to be preceeded by one space!
            else
                echo " GTR+I+G" >> $1/$2  # If error in modeltesting or parsing, set GTR+I+G as default model
                if [ $3 -eq 1 ]; then
                    echo -ne " WARN  | $(get_current_time) | No model extracted from file $file; setting model GTR+I+G as replacement.\n"
                fi
            fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No result files of the modeltesting (*.bestModel) found.\n"
        exit 1
    fi
}

modeltesting_via_sms()
#   This function conducts modeltesting via SMS for a series of input files; these input files start with the name "partition_"
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of SMS binary ($BinMDLTST)
#         $3: name of config file ($reformCFG)
#         $4: boolean variable if verbose or not ($verboseBOOL)
#   OUP:  none; generates *.bestModel* output files in temp folder
{
    # INTERNAL FUNCTION CHECKS
    (($# == 4)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done

    # USING CONFIG TEMPLATE
    CMDline=$(grep -A1 --ignore-case '^%SMS:' $1/$3 | tail -n1)
    CMDline=${CMDline//\[PATH_TO_SMS_SHELLSCRIPT\]/$2} # Substitution with Bash's built-in parameter expansion
    CMDline=$(eval echo $CMDline)
    # MODELTESTING
    count=$(ls -1 ./$1/partition_*.phy 2>/dev/null | wc -l)
    if [ $count != 0 ]; then
        for partit in ./$1/partition_*.phy; do
        partit_RUNTIME_start=$(date +%s)
        coreCMD=${CMDline//\[PARTITION_NAME\]/$partit} # Substitution with Bash's built-in parameter expansion; # NOTE: Don't overwrite the variable $CMDline to perform different substitutions in each loop, instead save it as a new variable ($coreCMD)
        coreCMD=${coreCMD//\[OUTPUT_FILE\]/${partit%.phy*}_sms.csv} # Substitution with Bash's built-in parameter expansion
        eval "$coreCMD 1>${partit%.phy*}.runlog 2>&1" # Executin of main_cmd
        partit_RUNTIME_dur=$(bc -l <<< "($(date +%s)-$partit_RUNTIME_start)/60")
        LC_ALL=C partit_RUNTIME_dur_pp=$(printf "%.3f minutes\n" $partit_RUNTIME_dur)
        if [ $4 -eq 1 ]; then
            echo -ne " INFO  | $(get_current_time) |   Processing complete for partition: $(basename $partit); Execution time: $partit_RUNTIME_dur_pp\n"
        fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No partition files (partition_*) found.\n"
        exit 1
    fi
}

########################################################################
########################################################################

# FUNCTIONS SPECIFIC TO 'PARTITIONFINDER'

assemble_mrbayes_block_from_partitionfinder()
#   This function ...
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of outfile ($NamOTFILE)
#   OUP:  append to the outfile
{
    # INTERNAL FUNCTION CHECKS
    (($# == 2)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    pf_bestscheme=$1/analysis/best_scheme.txt
    kw1="begin mrbayes;"
    kw2="end;"
    mrbayes_block=$(sed -n "/$kw1/,/$kw2/p" $pf_bestscheme)
    if [ ! -z "$mrbayes_block" ]; then
        echo "$mrbayes_block" >> $2 # NOTE: Append only! Also: $2 (i.e., NamOTFILE) does not reside in $1 (i.e., $tmpFOLDER)
    else
        echo -ne " ERROR | $(get_current_time) | Parsing of MRBAYES-block unsuccessful from file: $pf_bestscheme\n"
        exit 1
    fi
}

convert_datablock_to_phylip()
#   This function converts a NEXUS-formatted DATA-block to a PHYLIP-file.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: filename of DATA-block ($dataBlock)
#         $3: number of taxa in DNA matrix ($ntaxVar)
#         $4: number of characters in DNA matrix ($ncharVar)
#         $5: name of outfile ($phylipFil)
#   OUP:  writes a phylip-formatted file to disk
{
    # INTERNAL FUNCTION CHECKS
    (($# == 5)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # Write ntax and nvar to first line of phylip file
    echo "$3 $4" > $1/$5
    # Extract matrix from DATA-block
    pureMatrx=$(sed -n '/matrix/{:a;n;/;/b;p;ba}' $1/$2)
    #algnMatrx=$(echo "$pureMatrx" | column -t)
    #echo "$algnMatrx" >> $1/$5
    echo "$pureMatrx" >> $1/$5
}

extract_modelinfo_from_partitionfinder()
#   This function ...
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of file to collect bestfit model per partition ($bestModls)
#   OUP:  writes the bestModls file to disk
{
    # INTERNAL FUNCTION CHECKS
    (($# == 2)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    pf_bestscheme=$1/analysis/best_scheme.txt
    if [ -f "$pf_bestscheme" ]; then
        kw1="Subset \| Best Model"
        kw2="Scheme Description in PartitionFinder format"
        bestMdlOvw=$(sed -n "/$kw1/,/$kw2/p" $pf_bestscheme | grep '^[0-9]')
        bestModels=$(echo "$bestMdlOvw" | sed 's/[ \t]*//g' | awk -F'|' '{print $5" "$2}')
        if [ ! -z "$bestModels" ]; then
            echo "$bestModels" > $1/$2 # It's fine to overwrite the file here, as variable $bestModels already contains model names
        else
            echo -ne " WARN  | $(get_current_time) | Parsing of models unsuccessful from file: $pf_bestscheme\n"
        fi
    else
        echo -ne " ERROR | $(get_current_time) | No result files of the modeltesting (*.bestModel) found.\n"
        exit 1
    fi
}

modeltesting_via_partitionfinder()
#   This function conducts modeltesting via PartitionFinder
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: path to PartitionFinder Python script ($BinMDLTST)
#   OUP:  none; generates an output folder named 'analysis' in temp folder
{
    # INTERNAL FUNCTION CHECKS
    (($# == 2)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    # Test if Python2.7 present
    if command -v python2.7 1>/dev/null 2>&1; then
        :
    else
        echo -ne " ERROR | $(get_current_time) | Cannot find Python2.7 interpreter.\n"
        exit 1
    fi
    
    # Conduct modeltesting
    if [ -f "$1/partition_finder.cfg" ]; then
        python2 $2 $1 1>$1/partitionfinder_run.log 2>&1
    else
        echo -ne " ERROR | $(get_current_time) | No PartitionFinder config file (partition_finder.cfg) found.\n"
        exit 1
    fi
}

write_partitionfinder_config_file()
#   This function generates a config file as required by PartitionFinder.
#   INP:  $1: path to temp folder ($tmpFOLDER)
#         $2: name of config file ($reformCFG)
#         $3: name of phylip-formatted file ($phylipFil)
#         $4: name of file with charsets lines ($charSetsL)
#         $5: boolean variable if userscheme only or not ($usrschmBOOL)
#   OUP:  writes a phylip-formatted file to disk
{
    # INTERNAL FUNCTION CHECKS
    (($# == 5)) || { printf " ERROR | $(get_current_time) | The following function received an incorrect number of arguments: ${FUNCNAME[0]}\n"; exit 1; }
    for ((i=1; i<=$#; i++)); do
        argVal="${!i}"
        [[ -z "${argVal// }" ]] && { printf " ERROR | $(get_current_time) | The following function received an empty input argument: ${FUNCNAME[0]}\n"; exit 2; }
    done
    
    charsetLines=$(cat $1/$4)
    # USING CONFIG TEMPLATE
    CMDline=$(grep -A1 --ignore-case '^%PARTITIONFINDER:' $1/$2 | tail -n1)
    CMDline=${CMDline//\[NAME_OF_ALIGNMENT\]/$3} # Substitution with Bash's built-in parameter expansion
    CMDline=${CMDline//\[LIST_OF_DATA_BLOCKS\]/$charsetLines} # Substitution with Bash's built-in parameter expansion
    # Setting userscheme-only search in PartitionFinder or not.
    if [ $5 -eq 1 ]; then
        echo -ne " INFO  | $(get_current_time) | Estimating model fit of user-defined partition scheme only\n"
        CMDline=$(echo "$CMDline" | sed 's/\[schemes\].*/[schemes]/') # Remove everything after keyword "[schemes]"
        CMDline+="\nsearch=user;" # Appending string to variable
        charsetNums=$(echo "$charsetLines" | awk '{print "("$2")"}' | tr -d '\n') # NOTE: tr-command removes newlines from list
        CMDline+="\nuserdefined=$charsetNums;" # Appending string to variable
    fi
    echo -e "$CMDline" > $1/partition_finder.cfg
}

########################################################################
########################################################################

## EVALUATING SYSTEM, BASH SHELL AND SCRIPT EXECUTION
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Evaluating system, bash shell and script execution"\\n
fi

# Print system details to log
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) |   System info: $(uname -sr 2>/dev/null), proc: $(uname -p 2>/dev/null), arch: $(uname -m 2>/dev/null), numcores: $(get_num_of_avail_cores)"\\n
    echo -ne " INFO  | $(get_current_time) |   Bash info: $(bash --version | head -n1 2>/dev/null)"\\n
fi

# Print all bash arguments used to start this script
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) |   Command used: $THIS_SCRIPT $(echo $@)"\\n
fi

########################################################################

## CHECKING INFILES
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Checking infiles\n"
fi

# Renaming input variables
nexINFILE=$OPT_A
cfgINFILE=$OPT_B
typMDLTST=$OPT_C
BinMDLTST=$OPT_D
NamOTFILE=$OPT_E

if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) |   Type of Modeltesting selected: $typMDLTST\n"
fi

# Define outfile namestem
baseName=$(basename $nexINFILE) # Using basename to strip off path
filenStem=${baseName%.nex*} # Using parameter expansion to remove file extension

# Define outfile names
reformNEX=${filenStem}_ReformNEXFile
reformCFG=${filenStem}_ReformCFGFile
dataBlock=${filenStem}_NexusDatBlock
setsBlock=${filenStem}_NexusSetBlock
unspltMtx=${filenStem}_UnsplitMatrix
bestModls=${filenStem}_BestFitModels
charSetsL=${filenStem}_CharsetsLines
lsetDefns=${filenStem}_LsetDefinitns
phylipFil=${filenStem}_PhylipFormatd # Bash string extension cannot handle dots (i.e., must be "myfile_phy" instead of "myfile.phy")

# Checking input and output file availability
check_inp_outp_availability $nexINFILE $cfgINFILE $BinMDLTST $NamOTFILE $get_current_time

# Make temporary folder
#tmpFOLDER=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
if [ "$typMDLTST" = "partitionfinder" ] &&  [ "$usrschmBOOL" -eq 1 ]; then
    tmpFOLDER=${filenStem}_${typMDLTST}_userscheme_${RUNTIME_start_pp}_runFiles
else
    tmpFOLDER=${filenStem}_${typMDLTST}_${RUNTIME_start_pp}_runFiles
fi
mkdir -p $tmpFOLDER

########################################################################

## RE-FORMATTING INPUT NEXUS FILE
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Re-formatting input NEXUS file\n"
fi

reformat_nexus_file $nexINFILE $tmpFOLDER/$reformNEX
reformat_config_file $cfgINFILE $typMDLTST $tmpFOLDER/$reformCFG

# Extract charsets
cat $tmpFOLDER/$reformNEX | grep 'charset' > $tmpFOLDER/$charSetsL
# Extract number of characters in DNA matrix
ncharVar=$(grep ^dimensions $tmpFOLDER/$reformNEX | sed -n 's/.*nchar=\([^;]*\).*/\1/p') # NOTE: Terminating the search for ntax=xx with a semi-colon means that nchar has to come as last attribute (i.e., after ntax) in dimensions line.
# Extract number of taxa
ntaxVar=$(grep ^dimensions $tmpFOLDER/$reformNEX | sed -n 's/.*ntax=\([^ ]*\).*/\1/p') # NOTE: Terminating the search for ntax=xx with a space means that ntax has to come before nchar in dimensions line.

########################################################################

## SPLITTING NEXUS FILE INTO BLOCKS
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Splitting NEXUS file into blocks\n"
fi

split_nexus_into_blocks $tmpFOLDER $reformNEX $dataBlock $setsBlock

########################################################################

## CONFIRMING THAT PARTITIONS DO NOT OVERLAP
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Confirming that partitions do not overlap\n"
fi

test_if_partitions_overlap $tmpFOLDER $setsBlock

########################################################################

## ENSURING THAT PARTITIONS FORM A CONTINUOUS RANGE
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Ensuring that partitions form a continuous range\n"
fi

ensure_partitions_form_continuous_range $tmpFOLDER $setsBlock $ncharVar

########################################################################

## SPLITTING MATRIX INTO PARTITIONS
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Splitting matrix into partitions\n"
fi

split_matrix_into_partitions $tmpFOLDER $dataBlock $setsBlock

########################################################################

## PREPARING FILES FOR MODELTESTING
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Preparing files for modeltesting\n"
fi

if [ "$typMDLTST" = "sms" ]; then
    for partit in ./$tmpFOLDER/partition_*; do
        partit_ncharVar=$(grep ^dimensions $partit | sed -n 's/.*nchar=\([^;]*\).*/\1/p') # Extract number of characters in partition
        partit_name=$(basename $partit) # Because convert_datablock_to_phylip() concatenates path and filename, I need to extract basename of infile here
        convert_datablock_to_phylip $tmpFOLDER $partit_name $ntaxVar $partit_ncharVar ${partit_name}.phy
    done
fi

if [ "$typMDLTST" = "partitionfinder" ]; then
    convert_datablock_to_phylip $tmpFOLDER $dataBlock $ntaxVar $ncharVar $phylipFil
    write_partitionfinder_config_file $tmpFOLDER $reformCFG $phylipFil $charSetsL $usrschmBOOL
fi

########################################################################

## CONDUCTING MODELTESTING
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Conducting modeltesting via: $typMDLTST\n"
fi

if [ "$typMDLTST" = "jmodeltest" ]; then
    modeltesting_via_jmodeltest $tmpFOLDER $BinMDLTST $reformCFG $verboseBOOL
elif [ "$typMDLTST" = "sms" ]; then
    modeltesting_via_sms $tmpFOLDER $BinMDLTST $reformCFG $verboseBOOL
elif [ "$typMDLTST" = "partitionfinder" ]; then
    modeltesting_via_partitionfinder $tmpFOLDER $BinMDLTST
fi

########################################################################

## EXTRACTING INFORMATION FROM MODELTESTING RESULTS
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Extracting information from modeltesting results\n"
fi

if [ "$typMDLTST" = "jmodeltest" ]; then
    extract_modelinfo_from_jmodeltest $tmpFOLDER $bestModls $verboseBOOL
    convert_models_into_lset $tmpFOLDER $bestModls $lsetDefns $charSetsL
elif [ "$typMDLTST" = "sms" ]; then
    extract_modelinfo_from_sms $tmpFOLDER $bestModls $verboseBOOL
    convert_models_into_lset $tmpFOLDER $bestModls $lsetDefns $charSetsL
elif [ "$typMDLTST" = "partitionfinder" ]; then
    extract_modelinfo_from_partitionfinder $tmpFOLDER $bestModls
fi

########################################################################

## ASSEMBLING OUTPUT FILE
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Assembling output file\n"
fi

echo -e "#NEXUS\n" > $NamOTFILE

# Reconstitute DATA-block as interleaved
reconstitute_datablock_as_interleaved $tmpFOLDER $dataBlock $NamOTFILE 

# Append the SETS-block, which is commented out
echo -e "\n[\n$(cat $tmpFOLDER/$setsBlock)\n]\n" >> $NamOTFILE

# Append info on best-fitting models
echo -e "\n[\nBest-fitting models identified:\n$(cat $tmpFOLDER/$bestModls)\n]\n" >> $NamOTFILE

if [ "$typMDLTST" != "partitionfinder" ]; then # jmodeltest and sms both require function assemble_mrbayes_block_general()
    assemble_mrbayes_block_general $tmpFOLDER $NamOTFILE $lsetDefns $charSetsL $reformCFG
elif [ "$typMDLTST" = "partitionfinder" ]; then
    assemble_mrbayes_block_from_partitionfinder $tmpFOLDER $NamOTFILE
fi

########################################################################

## CLEANING UP THE WORK DIRECTORY
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Cleaning up\n"
fi

if [ $keepflsBOOL -eq 0 ]; then
    rm -r $tmpFOLDER
fi

# Stop timing execution time
RUNTIME_dur=$(bc -l <<< "($(date +%s)-$RUNTIME_start)/60")
LC_ALL=C RUNTIME_dur_pp=$(printf "%.3f minutes\n" $RUNTIME_dur)  # NOTE: "LC_ALL=C" is important due to different locales (German comma vs. world's decimal point)

# End of file message
if [ $verboseBOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Processing complete for file: $nexINFILE; Execution time: $RUNTIME_dur_pp\n"
fi

########################################################################

## EXIT WITHOUT MISTAKES
exit 0

########################################################################
