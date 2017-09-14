#!/bin/bash

########################################################################

THIS_SCRIPT=`basename ${BASH_SOURCE[0]}`
#DESCRIPTION="A bash shell scrip to convert a partitioned NEXUS file into a partitioned MrBayes analysis file"
AUTHOR="Michael Gruenstaeudl, PhD"
#COPYRIGHT="Copyright (C) 2015-2017 $AUTHOR"
#CONTACT="m.gruenstaeudl@fu-berlin.de"
VERSION="2017.09.14.1600"
USAGE="bash $THIS_SCRIPT -f <path_to_NEXUS_file> -m <path_to_modeltesting_tool> -t <type_of_modeltesting_tool>"

########################################################################

#    The script converts a DNA alignment in NEXUS format that contains character set ("charset") definitions into a partitioned NEXUS file ready for analysis with MrBayes. The character sets (hereafter "partitions") are hereby extracted, passed to modeltesting with one of two different modeltesting software tools (jModelTest or PartitionFinder2) to identify the best-fitting model of nucleotide substitition and then concatenated again. Where necessary, the best-fitting nucleotide substitution models are converted to the specifications read by MrBayes. A command block for MrBayes is appended to the NEXUS file that integrates these partition definitions.

#    Several tests regarding the character set definitions are performed during the execution of the script. First, the script evaluates if the entire alignment is covered by partitions. If not, additional partitions are defined so that all partitions together form a continuous range from {1} to {total size of matrix}. Second, the script evaluates if any of the character definitions overlap with one another. If an overlap is detected, the script exists with an error. Consequently, the initial character set definitions of the NEXUS file do not need to cover the entire alignment length, but must not be overlapping.

#    ARGS:
#        NEXUS file:    file path to, and name of, the NEXUS file
#        MODELTESTING tool: file path to, and name of, the modeltesting software tool
#        modeltesting TYPE: name of modeltesting software tool (available: jmodeltest, partitionfinder2)
#    OUTP:
#        MRBAYES file:  a NEXUS file ready for analysis with MrBayes

########################################################################

# PLANNED IMPROVEMENTS IN FUTURE VERSIONS


########################################################################

# NEW FEATURES:

# (a) Output now displays interleaved matrix.
# ORIGINAL PROBLEM: The Datablock in the output must not be copied from the input, but must be assembled as interleaved matrices from the individual partitions. The reason for this: MrBayes only accepts lines with 99990 characters, which may be smaller than the complete matrix. Thus, an interleaving is necessary.

#> It looks like the maximum allowed length of a token is defined in the 
#> source file mb.h.  You might try changing the value in the line:
#>
#> #define    CMD_STRING_LENGTH        100000
#>
#> to some value higher than the number of characters per line in your 
#> data set, then recompiling.

# (b) If an error during modeltesting for a single partition occurs (or the model is not correctly parsed), the code simply assigns a general model (GTR+I+G) to said partition.
# ORIGINAL PROBLEM: If an error during modeltesting for a single partition occurs, it should not crash the execution of the script for all partitions. Hence, future improvements of the code should make sure that an error during one modeltesting iterative (i.e., for one of several partitions) simply assigns a general model (GTR+I+G) to said partition if an error is thrown (i.e., comparable to a try-catch statement).

########################################################################

# INITIAL SETTINGS

# Toggle on: throw errors
set -e

# Start timing execution time
totalrun_time_start=$(date +%s)
totalrun_time_start_prettyp=$(date '+%Y.%m.%d.%H%M')

# Initialize variables to default values
OPT_A="file path to, and name of, the NEXUS file"
OPT_B="file path to, and name of, the jModelTest jar file"
keepBool=0
vrbseBool=0

# Set fonts for help function
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# Help function
function my_help {
    echo -e \\n"Help documentation for ${BOLD}$THIS_SCRIPT${NORM}."\\n
    echo -e "Version: $VERSION | Author: $AUTHOR"\\n
    echo -e "${REV}Usage:${NORM} $USAGE"\\n
    echo "Mandatory command line switches:"
    echo "${REV}-f${NORM}  --Sets file path to, and name of, ${BOLD}NEXUS file${NORM}. No default exists."
    echo "${REV}-m${NORM}  --Sets file path to, and name of, ${BOLD}modeltesting software tool${NORM}. No default exists."
    echo "${REV}-t${NORM}  --Sets ${BOLD}type${NORM} of modeltesting software tool (available: jmodeltest, partitionfinder2). No default exists."
    echo "Optional command line switches:"
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

while getopts :f:m:t:hkv FLAG; do
  case $FLAG in
    f)  OPT_A=$OPTARG
        ;;
    m)  OPT_B=$OPTARG
        ;;
    t)  OPT_C=$OPTARG
        ;;
    h)  my_help
        ;;
    k)  keepBool=1
        ;;
    v)  vrbseBool=1
        ;;
    \?) echo 'Invalid option: -$OPTARG' >&2
        exit 1
        ;;
    :)  echo 'Option -$OPTARG requires an argument.' >&2
        exit 1
        ;;
  esac
done

shift $((OPTIND-1))

########################################################################
########################################################################

# FUNCTIONS

get_current_time() {
    date '+%Y-%m-%d %Z %H:%M:%S'
}

reformat_nexus_file()
#   This function performs the following steps:
#   (a) it removes comments from a NEXUS file,
#   (b) it removes blank lines (both blank by whitespaces and by tabs), 
#   (c) it standardizes critical command lines (i.e. begin data, 
#       begin sets, dimensions, etc.) as lowercase, and
#   (d) it removes lines from the sets-block that do not start with the keyword "charset"
#   INP:  $1: name of complete NEXUS file
#         $2: name of a reformatted NEXUS file
#   OUP:  the re-formatted NEXUS file
{
    # Removes comments from a NEXUS file (i.e., delete all info enclosed by square brackets)
    nexWoCmts=$(sed 's/\[[^]]*\]//g' $1)
    # Removes blank lines (both blank by whitespaces and by tabs)
    nexWoBlkL=$(echo "$nexWoCmts" | sed '/^\s*$/d')
    # Standardizes critical command lines as lowercase (i.e., converts 
    # ENTIRE LINE that starts with a keyword TO LOWERCASE)
    # NOTE: These conversions to lowercase are critical for the correct 
    #       section identifikation below.
    # keyword1: "begin", line ends with semi-colon
    # NOTE: the following command converts "Begin Data;" to "begin data;", not just "begin"!
    reformKw1=$(echo "$nexWoBlkL" | awk 'BEGIN{IGNORECASE=1} /^ *begin\>.*; *$/ {$0=tolower($0)} 1')
    # keyword2: "end;"
    reformKw2=$(echo "$reformKw1" | awk 'BEGIN{IGNORECASE=1} /^ *end\;/ {$0=tolower($0)} 1')
    # keyword3: "matrix"
    reformKw3=$(echo "$reformKw2" | awk 'BEGIN{IGNORECASE=1} /^ *matrix\>/ {$0=tolower($0)} 1')
    # keyword4: "dimensions", line ends with semi-colon
    reformKw4=$(echo "$reformKw3" | awk 'BEGIN{IGNORECASE=1} /^ *dimensions\>.*; *$/ {$0=tolower($0)} 1')
    # keyword5: "format", line ends with semi-colon
    reformKw5=$(echo "$reformKw4" | awk 'BEGIN{IGNORECASE=1} /^ *format\>.*; *$/ {$0=tolower($0)} 1')
    # Convert only SPECIFIC KEYWORD TO LOWERCASE
    reformKw6=$(echo "$reformKw5" | awk 'tolower($1)=="charset"{$1=tolower($1)}1')
    # Removes lines from the sets-block that do not start with the keyword "charset"
    # NOTE: This check must come after the conversion to lowercase chars!
    # NOTE: This step generates the reformatted NEXUS file
    echo -e "#NEXUS\n" > $2
    echo "$reformKw6" | sed -n '/begin data\;/,/end\;/p' >> $2
    echo -e "\nbegin sets;" >> $2
    echo "$reformKw6" | sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' | grep 'charset' >> $2
    echo "end;" >> $2
    # Check that datatype is "DNA|dna"
    # NOTE: This check must come after the conversion to lowercase chars!
    dataTypeB=$(cat $2 | grep format | sed -n 's/.*datatype=\([^ ]*\).*/\1/p') # Extract info on datatype
    if [ ! "$dataTypeB" = "dna" ]; then
        echo -ne " ERROR | $(get_current_time) | Datatype not set as: DNA\n"
        exit 1
    fi
}

convert_models_into_lset()
#   This function converts the names of nucleotide substitution models 
#   into lset definitions readable by MrBayes.
#   NOTE: It also adds the model names as comments to the end of each lset line to allow the applyto-assignment (at the end of the function)!
#   INP:  $1: path to temp folder
#         $2: name of best-models file
#         $3: name of lset definitions file
#         $4: list of charsets
#   OUP:  none; operates on lset definitions file
#   NOTE:   The model abbreviation (e.g., `GTR+I+G`) must be at the line 
#           end and no whitespace must come thereafter. This end of line 
#           position precludes that string `GTR+I+G` would be replaced by 
#           `GTR+I`.
{
    # Step 1: Copy best-models file as lset definitions file while adding a comment sign to each line
    cp $1/$2 $1/$3
    awk -F'/' '{print "[# "$NF}' $1/$2 > $1/$3
    
    # Step 2: Edit the lset definitions file
    # NOTE: The final ampersand copies the model name to the end of the line; this is important to allow the applyto-assignment!
    ## GTR
    sed -i 's/.* GTR+I+G$/lset applyto=() nst\=6 rates\=invgamma\; &/' $1/$3
    sed -i 's/.* GTR+G$/lset applyto=() nst\=6 rates\=gamma\; &/' $1/$3
    sed -i 's/.* GTR+I$/lset applyto=() nst\=6 rates\=propinv\; &/' $1/$3
    sed -i 's/.* GTR$/lset applyto=() nst\=6\; &/' $1/$3
    ## SYM
    awk '/SYM\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* SYM+I+G$/lset applyto=() nst\=6 rates\=invgamma\; &/' $1/$3 # `/prset/!` means conduct replacement unless keyword ´prset´ present in line
    awk '/SYM\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* SYM+G$/lset applyto=() nst\=6 rates\=gamma\; &/' $1/$3
    awk '/SYM\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* SYM+I$/lset applyto=() nst\=6 rates\=propinv\; &/' $1/$3
    awk '/SYM$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* SYM$/lset applyto=() nst\=6\; &/' $1/$3
    ## HKY
    sed -i 's/.* HKY+I+G$/lset applyto=() nst\=2 rates\=invgamma\; &/' $1/$3
    sed -i 's/.* HKY+G$/lset applyto=() nst\=2 rates\=gamma\; &/' $1/$3
    sed -i 's/.* HKY+I$/lset applyto=() nst\=2 rates\=propinv\; &/' $1/$3
    sed -i 's/.* HKY$/lset applyto=() nst\=2\; &/' $1/$3
    ## K80(=K2P)
    awk '/K80\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* K80+I+G$/lset applyto=() nst\=2 rates\=invgamma\; &/' $1/$3
    awk '/K80\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* K80+G$/lset applyto=() nst\=2 rates\=gamma\; &/' $1/$3
    awk '/K80\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* K80+I$/lset applyto=() nst\=2 rates\=propinv\; &/' $1/$3
    awk '/K80$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* K80$/lset applyto=() nst\=2\; &/' $1/$3
    ## F81
    sed -i 's/.* F81+I+G$/lset applyto=() nst\=1 rates\=invgamma\; &/' $1/$3
    sed -i 's/.* F81+G$/lset applyto=() nst\=1 rates\=gamma\; &/' $1/$3
    sed -i 's/.* F81+I$/lset applyto=() nst\=1 rates\=propinv\; &/' $1/$3
    sed -i 's/.* F81$/lset applyto=() nst\=1\; &/' $1/$3
    ## JC
    awk '/JC\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* JC+I+G$/lset applyto=() nst\=1 rates\=invgamma\; &/' $1/$3
    awk '/JC\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* JC+G$/lset applyto=() nst\=1 rates\=gamma\; &/' $1/$3
    awk '/JC\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* JC+I$/lset applyto=() nst\=1 rates\=propinv\; &/' $1/$3
    awk '/JC$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3
    sed -i '/prset/! s/.* JC$/lset applyto=() nst\=1\; &/' $1/$3
    # Add a closing square bracket to all lines with `[#`
    awk '/\[\#/{$0 = $0 "]"} 1' $1/$3 > $1/tmp && mv $1/tmp $1/$3

    # Step 3: Replacing the partition names in the lset definitions with the line numbers in the charsets
    # NOTE: THIS IS WHY ORDERED PARTITIONS ARE CRITICAL!
    for charset in $(echo "$4" | awk '{print $2}'); do # Use doublequote-enclosing around a variable to maintain its newlines!
        lineNumbr=$(echo "$4" | awk '/'$charset'/{print NR; exit}')
        sed -i "/$charset/ s/applyto\=()/applyto\=($lineNumbr)/" $1/$3 # Double quotes are critical here due to variable
    done
}

ensure_partitions_form_continuous_range()
#   This function tests if the partitions form a continuous range. 
#   If the partitions don't form such a continuous range, additional
#   ranges are inserted such that all ranges taken together form a 
#   continuous range from {1} to {total size of matrix}.
#   INP:  $1: name of file containing the SETS-block
#         $2: total size of matrix
#   OUP:  update of $1
{
    # Get charset definition lines
    chrsetLns=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $1 | grep 'charset')
    # Convert from discrete to continuous range
    contRnge1=$(echo "$chrsetLns" | awk '{ split($4,curr,/[-;]/); currStart=curr[1]; currEnd=curr[2] } currStart > (prevEnd+1) { print "charset " "new"++cnt " = " prevEnd+1 "-" currStart-1 ";" } { print; prevEnd=currEnd }')
    # Add concluding partition, if missing
    matrxLngth=$2
    contRnge2=$(echo "$contRnge1" | tail -n1 | awk -F'[[:space:]]*|-' -v lngth=$matrxLngth ' $5<lngth {print "charset newFinal = " $5+1 "-" lngth ";"}')
    # Update SETS-block
    echo "begin sets;" > $1
    echo "$contRnge1" >> $1
    #echo -ne "$contRnge2" >> $1 # Option -ne necessary for pretty formatting only.
    echo "$contRnge2" >> $1
    echo "end;" >> $1
}

extract_modelinfo_from_jmodeltest()
#   This function extracts modelinfo from the output of jModeltest; the output was saved into logfiles ending with the name ".bestModel"
#   INP:  $1: path to temp folder ($tempFoldr)
#         $2: name of model overview file ($bestModls)
#         $3: boolean variable if verbose or not ($vrbseBool)
#         $4: function to get current time ($get_current_time)
#   OUP:  writes output to model overview file
{
    #count=`ls -1 $1/*.bestModel 2>/dev/null | wc -l`
    count=$(ls -1 $1/*.bestModel 2>/dev/null | wc -l)
    if [[ $count != 0 ]];
    then
        for file in ./$1/*.bestModel; do 
            echo -ne "$(basename $file)" | sed 's/partition_//g' | sed 's/\.bestModel//g' >> $1/$2
            bestModel=$(cat "$file" | grep -A1 ' Model selected:' | tail -n1 | sed -n 's/.*Model \= \([^ ]*\).*/\1/p' 2>/dev/null) # Parse best model from jModeltest output
            if [[ ! -z "$bestModel" ]];
            then
                echo " $bestModel" >> $1/$2  # NOTE: model has to be preceeded by one space!
            else
                echo " GTR+I+G" >> $1/$2  # If error in modeltesting or parsing, set GTR+I+G as standard model
                if [ $3 -eq 1 ]; then
                    echo -ne " WARNING | $(get_current_time) | No model extracted from file $file; setting model GTR+I+G as replacement.\n"
                fi
            fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No result files of the modeltesting (*.bestModel) found.\n"
        exit 1
    fi
}

modeltesting_via_jmodeltest()
#   This function conducts modeltesting via jmodeltest for a series of input files; these input files start with the name "partition_"
#   INP:  $1: path to temp folder ($tempFoldr)
#         $2: name of jModeltest binary ($mdltstBin)
#         $3: boolean variable if verbose or not ($vrbseBool)
#         $4: function to get current time ($get_current_time)
#   OUP:  none; generates *.bestModel* output files in temp folder
{
    #count=`ls -1 $1/partition_* 2>/dev/null | wc -l`
    count=$(ls -1 $1/partition_* 2>/dev/null | wc -l)
    if [[ $count != 0 ]];
    then
        for file in ./$1/partition_*; do
        partition_time_start=$(date +%s)
        # java -jar $2 -tr $nmbrCores -d "$file" -f -i -g 4 -AIC -o ${file}.bestModel 2>${file}.bestModel.err
        java -jar $2 -d "$file" -f -i -g 4 -AICc 1>${file}.bestModel 2>&1
        # JMODELTEST OPTIONS SELECTED:
        # -d sequenceFileName
        # -f include models with unequals base frecuencies
        # -i include models with a proportion invariable sites
        # -AIC calculate the Akaike Information Criterion
        # -g numberOfCategories: include models with rate variation among sites and number of categories
        # (-tr numberOfThreads)  # NOTE: Option no longer used, because: By default, the total number of cores in the machine is used by jModelTest (see p. 11 of jModelTest 2 Manual v0.1.10).
        partition_time_dur=$((($(date +%s)-$partition_time_start)/60))
        partition_time_dur_prettyp=$(printf "%.3f minutes\n" $partition_time_dur)
        if [ $3 -eq 1 ]; then
            echo -ne " INFO  | $(get_current_time) |   Processing complete for partition: $(basename $file); Execution time: $partition_time_dur_prettyp\n"
        fi
        done
    else
        echo -ne " ERROR | $(get_current_time) | No partition files (partition_*) found.\n"
        exit 1
    fi
}

reconstitute_datablock_as_interleaved()
#   This function reconstitutes a data block from individual partitions, but makes the data block interleaved.
#   INP:  $1: path to temp folder ($tempFoldr)
#         $2: name of DATA-block ($dataBlock)
#         $3: name of outfile ($outFilenm)
#   OUP:  none; writes to outfile
{
    # Step 1: Append dimension and format info of DATA-block
    cat $1/$2 | sed -n '/begin data\;/,/matrix/p' >> $3
    # Step 2: Specify that matrix is interleaved
    sed -i '/format/ s/\;/ interleave\=yes\;/' $3 # in-line replacement of semi-colon with interleave-info plus semi-colon, if line has keyword 'format'
    # Step 3: Append the individual partitions
    for p in $(ls $1/partition_* | grep -v bestModel); do
        echo -ne "\n[$(basename $p)]\n" >> $3
        pureMatrx=$(cat $p | sed -n '/matrix/{:a;n;/;/b;p;ba}')
        algnMatrx=$(echo "$pureMatrx" | column -t)
        echo "$algnMatrx" >> $3 # Append only the matrix of a partition, not the preceeding dimension and format info; also, don't append the closing ';\nend;' per partition
    done
    # Step 4: Append a closing ';\nend;'
    echo -e "\n;\nend;" >> $3
}

split_matrix_into_partitions()
#   This function takes the matrix of a NEXUS file and extracts from it a 
#   character range (hereafter "sub-matrix") specified by the definition 
#   of a charsets line.
#   INP:  $1: path to temp folder
#         $2: name of file containing the DATA-block
#         $3: name of file containing the SETS-block
#   OUP:  file with name $charrngFn in temp folder
#         file with name $partFname in temp folder
{
    pureMatrx=$(sed -n '/matrix/{:a;n;/;/b;p;ba}' $2)
    chrsetLns=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $3 | grep 'charset')
    nexusNew1='#NEXUS\n\n'
    nexusNew2=$(sed -e '/matrix/,$d' $2) # Get section between "#NEXUS" and "MATRIX"
    nexusNew3='\nmatrix'
    nexusNew4=';\nend;\n'
    myCounter=0
    while IFS= read -r line; do 
        myCounter=$((myCounter+1))
        partitFn1=$(printf %03d $myCounter)
        partitFn2=$(echo "$line" | awk '{print $2}')
        partitnFn=partition_${partitFn1}_${partitFn2}
        # Get the info on the charset range
        charstRng=$(echo "$line" | awk '{print $4}' | sed 's/\;//')
        # Step 1 of assembling the new partition file
        echo -e "$nexusNew1$nexusNew2$nexusNew3" > $1/$partitnFn
        # Step 2 of assembling the new partition file: Add the sub-matrix
        awk 'NR==FNR{start=$1;lgth=$2-$1+1;next} {print $1, substr($2,start,lgth)}' FS='-' <(echo "$charstRng") FS=' ' <(echo "$pureMatrx") >> $1/$partitnFn
        # Step 3 of assembling the new partition file
        echo -e "$nexusNew4" >> $1/$partitnFn
        # Get the length of the sub-matrix
        mtrxLngth=$(awk '{print $2-$1+1}' FS='-' <(echo "$charstRng"))
        # Replace the number of characters with the length of the sub-matrix
        sed -i "/dimensions / s/nchar\=.*\;/nchar\=$mtrxLngth\;/" $1/$partitnFn
    done <<< "$chrsetLns" # Using a here-string
}

split_nexus_into_blocks()
#   This function splits a NEXUS file into individual blocks.
#   INP:  $1: name of complete NEXUS file
#         $2: name of file containing the extracted DATA-block
#         $3: name of file containing the extracted SETS-block
#   OUP:  file with DATA-block
#         file with SETS-block
{
    # Extract the blocks
    cat $1 | sed -n '/begin data\;/,/end\;/p' > $2
    cat $1 | sed -n '/begin sets\;/,/end\;/p' > $3
}

test_if_partitions_overlap()
#   This function tests if any of the partitions (i.e., charset ranges)
#   overlap. If they do, an error is thrown.
#   INP:  $1: name of file containing the SETS-block
#   OUP:  none
{
    chrsetLns=$(sed -n '/begin sets\;/{:a;n;/end\;/b;p;ba}' $1 | grep 'charset')
    charstRng=$(echo "$chrsetLns" | awk '{print $4}' | sed 's/\;//')
    # Test if any of the charset ranges are overlapping
    charstOlp=$(echo "$charstRng" | awk -F"-" 'Q>=$1 && Q{print val}{Q=$NF;val=$0}')
    if [ ! "$charstOlp" = "" ]; then 
        echo -ne " ERROR | $(get_current_time) | Charset range overlaps with subsequent range: $charstOlp\n"
        exit 1
    fi
}

write_mrbayes_block()
#   This function appends a MrBayes block to a file ($outFilenm).
#   INP:  $1: path to temp folder
#         $2: name of outfile
#         $3: name of lset definitions file
#         $4: charsets variable
#   OUP:  none; operates on lset definitions file
{
    echo -e 'begin mrbayes;' >> $2
    echo -e 'set autoclose=yes;' >> $2
    echo -e 'set nowarnings=yes;' >> $2
    echo "$4" >> $2
    echo -n 'partition combined =' $(echo "$4" | wc -l) ': ' >> $2 # Line must not end with linebreak, thus -n
    echo "$4" | awk '{print $2}' | awk '{ORS=", "; print; }' >> $2
    sed -i 's/\(.*\)\,/\1;/' $2 # Replace final comma with semi-colon; don't make this replacement global
    echo -e '\nset partition = combined;' >> $2 # Option -e means that \n generates new line
    cat $1/$3 >> $2
    echo 'prset applyto=(all) ratepr=variable;' >> $2
    echo 'unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);' >> $2
    echo -e 'mcmcp ngen=20000000 temp=0.1 samplefreq=10000;\nmcmc;' >> $2
    echo -e 'end;\n\nquit;' >> $2
}

########################################################################
########################################################################

## STEP 00: CHECKING SYSTEM AND BASH SHELL
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 00: Evaluate system and bash shell\n"
fi

# Evaluate number of cores available
nmbrCores=$(nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null) # The "||" indicates that if "nproc" fails, do "grep ^processor /proc/cpuinfo".
posIntgr='^[0-9]+$'
if ! [[ $nmbrCores =~ $posIntgr ]]; then
    nmbrCores=1
fi

# Print system details to log
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) |   System info: $(uname -sr 2>/dev/null), proc: $(uname -p 2>/dev/null), arch: $(uname -m 2>/dev/null), numcores: $nmbrCores\n"
    echo -ne " INFO  | $(get_current_time) |   Bash info: $(bash --version | head -n1 2>/dev/null)\n"
fi

########################################################################

## STEP 01: CHECKING INFILES
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 01: Checking infiles\n"
fi

# Renaming input variables
nexusFile=$OPT_A
mdltstBin=$OPT_B
mdltstTyp=$OPT_C

# Checking if input files exists
if [[ ! -f $nexusFile ]]; then 
    echo -ne " ERROR | $(get_current_time) | Not found: $nexusFile\n"
    exit 1
fi
if [[ ! -f $mdltstBin ]]; then 
    echo -ne " ERROR | $(get_current_time) | Not found: $mdltstBin\n"
    exit 1
fi

# Define outfile namestem
baseName=$(basename $nexusFile) # Using basename to strip off path
filenStem=${baseName%.nex*} # Using parameter expansion to remove file extension

# Define outfile names
reformNex=${filenStem}_ReformInpFile
dataBlock=${filenStem}_NexusDatBlock
setsBlock=${filenStem}_NexusSetBlock
unspltMtx=${filenStem}_UnsplitMatrix
bestModls=${filenStem}_BestFitModels
lsetDefns=${filenStem}_LsetDefinitns
outFilenm=${filenStem}_${totalrun_time_start_prettyp}_partitioned.mrbayes

# Check if outfile already exists
if [[ $outFilenm = /* ]];
then
    echo -ne " ERROR | $(get_current_time) | Outfile already exists in filepath: $outFilenm\n"
    exit 1
fi

# Make and cd into temporary folder
tempFoldr=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
mkdir -p $tempFoldr

########################################################################

## STEP 02: RE-FORMATTING INPUT NEXUS FILE
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 02: Re-formatting input NEXUS file\n"
fi

reformat_nexus_file $nexusFile $tempFoldr/$reformNex

# Extract charsets
charSetsL=$(cat $tempFoldr/$reformNex | grep 'charset')

########################################################################

## STEP 03: SPLITTING NEXUS FILE INTO BLOCKS
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 03: Splitting NEXUS file into blocks\n"
fi

split_nexus_into_blocks $tempFoldr/$reformNex $tempFoldr/$dataBlock $tempFoldr/$setsBlock

# Check if dataBlock successfully generated
if [ ! -s $tempFoldr/$dataBlock ];
then
    echo -ne " ERROR | $(get_current_time) | No file size: $dataBlock\n"
    exit 1
fi

# Check if setsBlock successfully generated
if [ ! -s $tempFoldr/$setsBlock ];
then
    echo -ne " ERROR | $(get_current_time) | No file size: $setsBlock\n"
    exit 1
fi

########################################################################

## STEP 04: CONFIRM THAT PARTITIONS DO NOT OVERLAP
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 04: Confirming that partitions do not overlap\n"
fi

test_if_partitions_overlap $tempFoldr/$setsBlock

########################################################################

## STEP 05: ENSURE THAT PARTITIONS FORM A CONTINUOUS RANGE
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 05: Ensuring that partitions form a continuous range\n"
fi

# Extract number of characters in DNA matrix
ncharVar=$(grep ^dimensions $tempFoldr/$dataBlock | sed -n 's/.*nchar=\([^;]*\).*/\1/p')

ensure_partitions_form_continuous_range $tempFoldr/$setsBlock $ncharVar

########################################################################

## STEP 06: SPLITTING MATRIX INTO PARTITIONS
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 06: Splitting matrix into partitions\n"
fi

split_matrix_into_partitions $tempFoldr $tempFoldr/$dataBlock $tempFoldr/$setsBlock

########################################################################

## STEP 07: CONDUCTING MODELTESTING
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 07: Conducting modeltesting\n"
fi

modeltesting_via_jmodeltest $tempFoldr $mdltstBin $vrbseBool $get_current_time

########################################################################

## STEP 08: EXTRACTING MODEL INFORMATION FROM JMODELTEST2
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 08: Extracting model information from jModelTest2\n"
fi

extract_modelinfo_from_jmodeltest $tempFoldr $bestModls $vrbseBool $get_current_time

########################################################################

## STEP 09: CONVERTING MODELTEST RESULTS TO LSET SPECS FOR MRBAYES
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 09: Converting modeltest results to lset definitions for MrBayes\n"
fi

if [[ ! -f $tempFoldr/$bestModls ]];
then 
    echo -ne " ERROR | $(get_current_time) | File not found: $tempFoldr/$bestModls\n"
    exit 1
fi

# Setting up lset specifications
convert_models_into_lset $tempFoldr $bestModls $lsetDefns "$charSetsL"
# NOTE: $charSetsL is a multi-line variable, whose linebreaks shall be maintained; thus, it is passed as doublequote-enclosed.

########################################################################

## STEP 10: ASSEMBLING OUTPUT FILE
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 10: Assembling partitions into output file\n"
fi

echo -e "#NEXUS\n" > $outFilenm

# Reconstitute DATA-block as interleaved
reconstitute_datablock_as_interleaved $tempFoldr $dataBlock $outFilenm 

# Append the SETS-block, which is commented out
echo -e "\n[\n$(cat $tempFoldr/$setsBlock)\n]\n" >> $outFilenm

# Append info on best-fitting models
echo -e "\n[\nBest-fitting models identified:\n$(cat $tempFoldr/$bestModls)\n]\n" >> $outFilenm

########################################################################

## STEP 11: GENERATING MRBAYES BLOCK
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 11: Generating MrBayes block\n"
fi

write_mrbayes_block $tempFoldr $outFilenm $lsetDefns "$charSetsL"
# NOTE: $charSetsL is a multi-line variable, whose linebreaks shall be maintained; thus, it is passed as doublequote-enclosed.

########################################################################

## STEP 12: CLEANING UP THE WORK DIRECTORY
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 12: Cleaning up\\n"
fi

if [ $vrbseBool -eq 1 ]; then
    runFolder=${filenStem}_${totalrun_time_start_prettyp}_runFiles
    mkdir -p $runFolder
    cp -r $tempFoldr/* $runFolder
fi
rm -r $tempFoldr

# Stop timing execution time
totalrun_time_dur=$((($(date +%s)-$totalrun_time_start)/60))
totalrun_time_dur_prettyp=$(printf "%.3f minutes\n" $totalrun_time_dur)

# End of file message
if [ $vrbseBool -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Processing complete for file: $nexusFile; Execution time: $totalrun_time_dur_prettyp\n"
fi

########################################################################

# Exit without mistakes
exit 0

########################################################################
