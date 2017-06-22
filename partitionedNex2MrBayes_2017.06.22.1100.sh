#!/bin/bash

########################################################################

THIS_SCRIPT=`basename ${BASH_SOURCE[0]}`
#DESCRIPTION="A bash shell scrip to convert a partitioned NEXUS file into a partitioned MrBayes analysis file"
AUTHOR="Michael Gruenstaeudl, PhD"
#COPYRIGHT="Copyright (C) 2015-2017 $AUTHOR"
#CONTACT="m.gruenstaeudl@fu-berlin.de"
VERSION="2017.06.22.1100"
USAGE="bash $THIS_SCRIPT -f <path_to_NEXUS_file> -m <path_to_modeltest.jar>"

########################################################################

#    This script converts a DNA alignment in NEXUS format that contains 
#    character set (i.e., charset) specifications into a partitioned 
#    NEXUS file ready for analysis with MrBayes. The character sets are 
#    hereby extracted, passed to modeltesting with jModelTest to identify 
#    the best-fitting model of nucleotide substitition and then 
#    concatenated. Only regions of the input alignment with a character 
#    set designations are passed to modeltesting and are concatenated
#    into the output alignment. The best-fitting nucleotide subtsitution 
#    models are converted to the specifications read by MrBayes.    
#    
#    ARGS:
#        NEXUS file:    file path to, and name of, the NEXUS file
#        MODELTEST jar: file path to, and name of, the jModelTest jar file
#    OUTP:
#        MRBAYES file:  a NEXUS file ready for analysis with MrBayes

########################################################################

# TO-DO LIST

# (a)   Include a test to see if the entire alignment is covered by 
#       charsets. If not, spit out an ERROR.
# (b)   Include a test to see if the character definitions overlap. 
#       If not, spit out an ERROR. Charset definitions must not 
#       overlap. Otherwise regions will be counted twice.
# (c)   Design such that STEPS 3, 4 and 5 are integrated into the loop.
# (d)   Wrap the N-of-cores evaluation into a try-statement 
#       (so that N=1 can be selected in the worst scenario).
# (e)   Make the selection of the number of cores an input variable. 
#       The number of cores are only inferred via 
#       unless ´grep -c ^processor /proc/cpuinfo´ if the user has set 
#       them manually.

########################################################################

# INITIAL SETTINGS

# Exit on error
set -e

# Initialize variables to default values
OPT_A="file path to, and name of, the NEXUS file"
OPT_B="file path to, and name of, the jModelTest jar file"
VERBOSE_BOOL=0

# Set fonts for help function
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# Help function
function HELP {
    echo -e \\n"Help documentation for ${BOLD}$THIS_SCRIPT${NORM}."\\n
    echo -e "Version: $VERSION | Author: $AUTHOR"\\n
    echo -e "${REV}Usage:${NORM} $USAGE"\\n
    echo "Mandatory command line switches:"
    echo "${REV}-f${NORM}  --Sets file path to, and name of, ${BOLD}NEXUS file${NORM}. No default exists."
    echo "${REV}-m${NORM}  --Sets file path to, and name of, ${BOLD}jModelTest jar file${NORM}. No default exists."
    echo "Optional command line switches:"
    echo "${REV}-v${NORM}  --Displays ${BOLD}verbose${NORM} status messages. Default is ${BOLD}off${NORM}."
    echo -e "${REV}-h${NORM}  --Displays this help message."\\n
    exit 1
}

# Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
#echo -e \\n"Number of arguments: $NUMARGS"
if [ $NUMARGS -eq 0 ]; then
    HELP
fi

########################################################################

# GETOPTS CODE

while getopts :f:m:hv FLAG; do
  case $FLAG in
    f)  OPT_A=$OPTARG
        ;;
    m)  OPT_B=$OPTARG
        ;;
    h)  HELP
        ;;
    v)  VERBOSE_BOOL=1
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

# FUNCTIONS

get_current_time() {
    date '+%Y-%m-%d %H:%M:%S'
}

split_nexus_into_blocks()
#   This function splits a NEXUS file into individual blocks.
#   INP:  $1: name of complete NEXUS file
#         $2: name of file containing the extracted DATA-block
#         $3: name of file containing the extracted SETS-block
#         $4: name of file containing the unsplit matrix
#   OUP:  file with DATA-block
#         file with SETS-block
{
    # Note: A conversion of the NEXUS file to lower case only is not possible, because model matching would not work
    cat $1 | sed -n '/begin data\;\|BEGIN DATA\;\|Begin Data\;/,/end\;\|END\;\|End\;/p' > $2
    cat $1 | sed -n '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/,/end\;\|END\;\|End\;/p' > $3
    sed -n '/matrix\|MATRIX\|Matrix/{:a;n;/;/b;p;ba}' $2 > $4
}

confirm_validity_of_partitions()
#   This function ... .
#   It eventually sets up the file containing the corrected, unsplit charsets
#   INP:  $1: name of ...
#         $2: name of file containing the SETS-block
#         $3: name of file containing the corrected, unsplit charsets
#   OUP:  file with ...
{
    #sed -n '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/{:a;n;/end\;\|END\;\|End\;/b;p;ba}' $2 | grep 'charset\|CHARSET\|CharSet' > $3
    sed -n '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/{:a;n;/end\;\|END\;\|End\;/b;p;ba}' $1 | grep 'charset\|CHARSET\|CharSet' > $2
}

split_matrix_into_partitions()
#   This function takes the matrix of a NEXUS file and extracts from it a 
#   character range (hereafter "sub-matrix") specified by the definition 
#   of a charsets line.
#   INP:  $1: path to temp folder
#         $2: name of complete NEXUS file
#         $3: name of file containing the unsplit matrix
#         $4: name of file containing the unsplit charsets
#   OUP:  file with name $CHARRNG_FNAME in temp folder
#         file with name $PARTITION_FNAME in temp folder
{
    NEXUS_P1='\n#NEXUS'
    NEXUS_P2=$(sed -e '/matrix\|MATRIX\|Matrix/,$d' $2)
    NEXUS_P3='\nMATRIX'
    NEXUS_P4=';\nEND;\n'
    COUNTER=0
    while IFS= read -r line; do 
        COUNTER=$((COUNTER+1))
        CHARSET_FNAME_P1=$(printf %03d $COUNTER)
        CHARSET_FNAME_P2=$(echo "$line" | awk '{print $2}')             # Get the gene name, which is the second word per line.
        CHARSET_FNAME=${CHARSET_FNAME_P1}_${CHARSET_FNAME_P2}           # Prepend number to the charset name
        CHARRNG_FNAME=${CHARSET_FNAME}_range
        PARTITION_FNAME=partition_${CHARSET_FNAME}                      # Define partition filename
        echo "$line" | awk '{print $4}' | sed 's/\;//' > $1/$CHARRNG_FNAME  # Get the info on the charset range
        echo -e "$NEXUS_P1 $NEXUS_P2 $NEXUS_P3" > $1/$PARTITION_FNAME             # Step 1 of assembling the new partition file
        awk 'NR==FNR{start=$1;lgth=$2-$1+1;next} {print $1, substr($2,start,lgth)}' FS='-' $1/$CHARRNG_FNAME FS=' ' $3 >> $1/$PARTITION_FNAME  # Step 2 of assembling the new partition file: Add the sub-matrix
        echo -e "$NEXUS_P4" >> $1/$PARTITION_FNAME                      # Step 3 of assembling the new partition file
        NEW_LENGTH=$(awk '{print $2-$1+1}' FS='-' $TMPFLD/$CHARRNG_FNAME)   # Get the length of the sub-matrix
        sed -i "/dimensions\|DIMENSIONS\|Dimensions/ s/NCHAR\=.*\;/NCHAR\=$NEW_LENGTH\;/" $1/$PARTITION_FNAME  # Replace the number of characters with the length of the sub-matrix
    done < $4
}

convert_models_into_lset()
#   This function converts the names of nucleotide substitution models 
#   into lset definitions readable by MrBayes.
#   INP:  $1: path to temp folder
#         $2: name of lset specs file
#   OUP:  none; operates on lset specs file
#   NOTE:   The model abbreviation (e.g., `GTR+I+G`) must be at the line 
#           end and no whitespace must come thereafter. This end of line 
#           position precludes that string `GTR+I+G` is also replaced by 
#           `GTR+I`.
{
##GTR
    sed -i 's/.* GTR+I+G$/lset applyto=() nst\=6 rates\=invgamma\; &/' $1/$2
    sed -i 's/.* GTR+G$/lset applyto=() nst\=6 rates\=gamma\; &/' $1/$2
    sed -i 's/.* GTR+I$/lset applyto=() nst\=6 rates\=propinv\; &/' $1/$2
    sed -i 's/.* GTR$/lset applyto=() nst\=6\; &/' $1/$2
##SYM
    awk '/SYM\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* SYM+I+G$/lset applyto=() nst\=6 rates\=invgamma\; &/' $1/$2 # `/prset/!` means conduct replacement unless keyword ´prset´ present in line
    awk '/SYM\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* SYM+G$/lset applyto=() nst\=6 rates\=gamma\; &/' $1/$2
    awk '/SYM\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* SYM+I$/lset applyto=() nst\=6 rates\=propinv\; &/' $1/$2
    awk '/SYM$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* SYM$/lset applyto=() nst\=6\; &/' $1/$2
##HKY
    sed -i 's/.* HKY+I+G$/lset applyto=() nst\=2 rates\=invgamma\; &/' $1/$2
    sed -i 's/.* HKY+G$/lset applyto=() nst\=2 rates\=gamma\; &/' $1/$2
    sed -i 's/.* HKY+I$/lset applyto=() nst\=2 rates\=propinv\; &/' $1/$2
    sed -i 's/.* HKY$/lset applyto=() nst\=2\; &/' $1/$2
##K80(=K2P)
    awk '/K80\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* K80+I+G$/lset applyto=() nst\=2 rates\=invgamma\; &/' $1/$2
    awk '/K80\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* K80+G$/lset applyto=() nst\=2 rates\=gamma\; &/' $1/$2
    awk '/K80\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* K80+I$/lset applyto=() nst\=2 rates\=propinv\; &/' $1/$2
    awk '/K80$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* K80$/lset applyto=() nst\=2\; &/' $1/$2
##F81
    sed -i 's/.* F81+I+G$/lset applyto=() nst\=1 rates\=invgamma\; &/' $1/$2
    sed -i 's/.* F81+G$/lset applyto=() nst\=1 rates\=gamma\; &/' $1/$2
    sed -i 's/.* F81+I$/lset applyto=() nst\=1 rates\=propinv\; &/' $1/$2
    sed -i 's/.* F81$/lset applyto=() nst\=1\; &/' $1/$2
##JC
    awk '/JC\+I\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* JC+I+G$/lset applyto=() nst\=1 rates\=invgamma\; &/' $1/$2
    awk '/JC\+G$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* JC+G$/lset applyto=() nst\=1 rates\=gamma\; &/' $1/$2
    awk '/JC\+I$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* JC+I$/lset applyto=() nst\=1 rates\=propinv\; &/' $1/$2
    awk '/JC$/{$0 = $0 "\nprset applyto=() statefreqpr=fixed(equal);" FS $0} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
    sed -i '/prset/! s/.* JC$/lset applyto=() nst\=1\; &/' $1/$2
# Add a closing square bracket to all lines with `[#`
    awk '/\[\#/{$0 = $0 "]"} 1' $1/$2 > $1/tmp && mv $1/tmp $1/$2
}

write_mrbayes_block()
#   This function appends a MrBayes block to a file ($OUTFNAME).
#   INP:  $1: path to temp folder
#         $2: name of outfile
#         $3: name of lset specs file
#         $4: charsets variable
#   OUP:  none; operates on lset specs file
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
    echo -e "sump burnin=1000 filen=$NEXUS nrun=2;" >> $2 # Double quotes are critical here due to variable
    echo -e "sumt burnin=1000 filen=$NEXUS nrun=2;" >> $2 # Double quotes are critical here due to variable
    echo -e 'end;\n\nquit;' >> $2
}

########################################################################

## Step 01: Checking infiles
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 01: Checking infiles\n"
fi

# Renaming input variables
NEXUS=$OPT_A
MODELTEST=$OPT_B

# Checking if input files exists
if [[ ! -f $NEXUS ]]; then 
    echo -ne " ERROR | $(get_current_time) | Not found: $NEXUS\n"
    exit 1
fi
if [[ ! -f $MODELTEST ]]; then 
    echo -ne " ERROR | $(get_current_time) | Not found: $MODELTEST\n"
    exit 1
fi

# Define outfile namestem
BASENAME=$(basename $NEXUS) # Using basename to strip off path
FILESTEM=${BASENAME%.nex*} # Using parameter expansion to remove file extension

# Define outfile names
DATABLOC=${FILESTEM}_DataBlock
SETSBLOC=${FILESTEM}_SetsBlock
UNSPLT_M=${FILESTEM}_UnsplitMatrix
UNSPLT_C=${FILESTEM}_UnsplitCharsets
MDLOVRVW=${FILESTEM}_ModelOverview
LSETSPEC=${FILESTEM}_LsetSpecs
OUTFNAME=${FILESTEM}_partitioned.mrbayes

# Parsing charsets from NEXUS file
CHARSETS=$(cat $NEXUS | grep 'charset\|CharSet\|CHARSET')

# Evaluate number of cores available on Desktop machine
CORES=$(grep -c ^processor /proc/cpuinfo)

# Check if outfile already exists
if [[ $OUTFNAME = /* ]];
then
    echo -ne " ERROR | $(get_current_time) | Outfile already exists in filepath: $OUTFNAME\n"
    exit 1
fi

# Make and cd into temporary folder
TMPFLD=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
mkdir -p $TMPFLD

########################################################################

## Step 02: Splitting NEXUS file into blocks
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 02: Splitting NEXUS file into blocks\n"
fi

split_nexus_into_blocks $NEXUS $TMPFLD/$DATABLOC $TMPFLD/$SETSBLOC $TMPFLD/$UNSPLT_M

########################################################################

## Step 03a: Confirm validity of partitions
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 03a: Confirm validity of partitions\n"
fi

# CONTINUE HERE!
#exit 1
confirm_validity_of_partitions $TMPFLD/$SETSBLOC $TMPFLD/UNSPLT_C

########################################################################

## Step 03b: Splitting matrix into partitions
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 03b: Splitting matrix into partitions\n"
fi

split_matrix_into_partitions $TMPFLD $TMPFLD/$DATABLOC $TMPFLD/$UNSPLT_M $TMPFLD/UNSPLT_C

########################################################################

## Step 04: Conducting modeltesting via jModelTest2
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 04: Conducting modeltesting via jModelTest2\n"
fi

count=`ls -1 $TMPFLD/partition_* 2>/dev/null | wc -l`
if [[ $count != 0 ]];
then
    for file in ./$TMPFLD/partition_*; do 
    java -jar $MODELTEST -tr $CORES -d "$file" -g 4 -i -f -AIC -o ${file}.bestModel 1>${file}.bestModel.log 2>&1
    done
else
    echo -ne " ERROR | $(get_current_time) | No partition files (partition_*) found.\n"
    exit 1
fi

########################################################################

## Step 05: Extracting model information from jModelTest2
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 05: Extracting model information from jModelTest2\n"
fi

count=`ls -1 $TMPFLD/*.bestModel 2>/dev/null | wc -l`
if [[ $count != 0 ]];
then
    for file in ./$TMPFLD/*.bestModel; do 
    echo -ne "$file" | sed 's/partition_//g' | sed 's/\.bestModel//g' >> $TMPFLD/$MDLOVRVW
    cat "$file" | grep -A1 ' Model selected:' | tail -n1 | sed 's/Model = //g' | sed 's/   / /g' >> $TMPFLD/$MDLOVRVW  # Have to use command ´cat´, cannot use ´echo´ here; not sure why.
    done
else
    echo -ne " ERROR | $(get_current_time) | No result files of the modeltesting (*.bestModel) found.\n"
    exit 1
fi

########################################################################

## Step 06: Converting modeltest results to lset specs for MrBayes
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 06: Converting modeltest results to lset specs for MrBayes\n"
fi

if [[ ! -f $TMPFLD/$MDLOVRVW ]];
then 
    echo -ne " ERROR | $(get_current_time) | File not found: $TMPFLD/$MDLOVRVW\n"
    exit 1
fi

awk -F'/' '{print "[# "$NF}' $TMPFLD/$MDLOVRVW > $TMPFLD/$LSETSPEC

# Retrieving the list of best fitting models at this point
BEST_MDLS=$(sed 's/\[\# //g' $TMPFLD/$LSETSPEC)

# Setting up lset specifications
convert_models_into_lset $TMPFLD $LSETSPEC

########################################################################

## Step 07: Replacing partition names with partition numbers (This is why ordered partitions are CRITICAL!)
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 07: Replacing partition names with partition numbers\n"
fi

# Replacing the partition name sin $LSETSPEC with the line numbers in the NEXUS file
for file in $(echo "$CHARSETS" | awk '{print $2}'); do # Use doublequote-enclosing around a variable to maintain its newlines!
LINENUM=$(echo "$CHARSETS" | awk '/'$file'/{print NR; exit}');
sed -i "/$file/ s/applyto\=()/applyto\=($LINENUM)/" $TMPFLD/$LSETSPEC # Double quotes are critical here due to variable
done

########################################################################

## Step 08: Assembling output file
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 08: Assembling output file\n"
fi

echo -e "#NEXUS\n" > $OUTFNAME
cat $TMPFLD/$DATABLOC >> $OUTFNAME
echo -e "\n[" >> $OUTFNAME # SETS-block is commented out
cat $TMPFLD/$SETSBLOC >> $OUTFNAME
echo -e "]" >> $OUTFNAME

# Comment out SETS block
#sed -i '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/ s/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/[\nbegin sets\;/g' $OUTFNAME
#echo -e ']\n' >> $OUTFNAME

echo -e "\n[\nBest-fitting models identified:\n$BEST_MDLS\n]\n" >> $OUTFNAME

########################################################################

## Step 09: Generating MrBayes block
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 09: Generating MrBayes block\n"
fi

write_mrbayes_block $TMPFLD $OUTFNAME $LSETSPEC "$CHARSETS"

########################################################################

## Step 10: Cleaning up the work directory
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Step 10: Cleaning up\\n"
fi

rm -r $TMPFLD

# End of file message
if [ $VERBOSE_BOOL -eq 1 ]; then
    echo -ne " INFO  | $(get_current_time) | Processing complete for: $NEXUS\n"
fi

# Exit without mistakes
exit 0

########################################################################
