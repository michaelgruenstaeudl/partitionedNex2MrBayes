#!/bin/bash
set -e # Exit on error
CURRENTTIME() { date '+%Y-%m-%d %H:%M:%S'; }

########################################################################

#Initialize variables to default values.
NEXUS=""
MODELTEST=""

## GETOPTS CODE

while getopts :f:m:hv FLAG; do
  case $FLAG in
    f)  NEXUS=$OPTARG
        ;;
    m)  MODELTEST=$OPTARG
        ;;
    h)  echo 'Help instructions here!' >&2
        ;;
    v)  echo 'Verbose mode on!' >&2
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

## LEGACYCODE - DELETE, IF NOT USED

## Checking number of arguments (NEXUS file, ModelTestJar)
#if [[ $# != 2 ]];
#then
#    echo -ne " ERROR | $(CURRENTTIME) | Incorrect number of input files specified.\n"
#    exit 1
#fi
#
## Assign infiles to variables
#NEXUS=$1
#MODELTEST=$2 # /home/michael_science/binaries/jModelTest_2.1.7/jModelTest.jar

########################################################################

TITLE="nex2partMrBayes.sh"
DESCRIPTION="Shell script to convert NEXUS alignment with charsets to partitioned MrBayes analysis file"
AUTHOR="Michael Gruenstaeudl, PhD"
COPYRIGHT="Copyright (C) 2015-2017 Michael Gruenstaeudl"
CONTACT="m.gruenstaeudl@fu-berlin.de"
VERSION="2017.06.21.1700"
USAGE="bash nex2partMrBayes.sh <abs_path_to_NEXUS_file> <abs_path_to_modeltest.jar>"
DEPENDENCIES="python2, biopython"

echo -ne "\n Title: $TITLE | Version: $VERSION | Author: $AUTHOR\n"
echo -ne " Usage: $USAGE\n"

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
#    Args:
#        NEXUS file:   absolute path to, and name of, the NEXUS file
#        MODELTEST jar: absolute path to, and name of, the jModelTest jar file
#    Outp:
#        MRBAYES NEXUS file: a NEXUS file ready for analysis with MrBayes
#    
#    
#    Example:
#        sh /home/michael_science/git/michaelgruenstaeudl_PlastidGenomics/09_Nex2PartitionedMrBayes/
#            
#            Nex2PartitionedMrBayes_2017.05.10.1700.sh /home/michael_science/Desktop/TESTFELD/partitioned_analysis/Nympheal_CompleteCDS_15Taxa_2017.02.08_CharsetsOrdered.nex
#            
#            /home/michael_science/binaries/jModelTest_2.1.7/jModelTest.jar
#    
#    
#        python2 ReassemblePlastome_2017.02.08.1800.py
#            -i TESTFELD/Nymphaea_ampla_KU189255_GeneNamesHomogenized.gb
#            -f /home/michael_science/Desktop/TESTFELD/05_collected/
#            -s 89952
#            -e 115128
#            -o TESTFELD/GeneNames.txt

########################################################################

## TO DO

# - Include a test to see if the entire alignment is covered by charsets. If not, spit out an ERROR. 
#
# - Include a test to see if the character definitions overlap. If not, spit out an ERROR. Charset definitions must not overlap. Otherwise regions will be counted twice.
#
# - Provide a -h help command that provides all detailed info. Also: --verbose, --version, --info
#
# - Wrap the N-of-cores evaluation into a try-statement (so that N=1 can be selected in the worst scenario).
#
# - Design such that STEPS 3, 4 and 5 are integrated into the loop.
#
# - Make the selection of the number of cores an input variable. The number of cores are only inferred via ´grep -c ^processor /proc/cpuinfo´ if the user hasn't set it manually.
#

########################################################################

## FUNCTIONS

split_nexus_into_blocks()
#   This function splits a NEXUS file into individual blocks and sub-blocks.
#   INP:  $1: complete NEXUS file
#   OUP:  file with DATA-block
#         file with SETS-block
#         file with DNA matrix only
#         file with charset lines only
{
    # Note: A conversion of the NEXUS file to lower case only is not possible, because model matching would not work
    cat $NEXUS | sed -n '/begin data\;\|BEGIN DATA\;\|Begin Data\;/,/end\;\|END\;\|End\;/p' > $TMPFLD/$DATABLOC
    cat $NEXUS | sed -n '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/,/end\;\|END\;\|End\;/p' > $TMPFLD/$SETSBLOC
    sed -n '/matrix\|MATRIX\|Matrix/{:a;n;/;/b;p;ba}' $TMPFLD/$DATABLOC > $TMPFLD/unsplit_matrix
    sed -n '/begin sets\;\|BEGIN SETS\;\|Begin Sets\;/{:a;n;/end\;\|END\;\|End\;/b;p;ba}' $TMPFLD/$SETSBLOC | grep 'charset\|CHARSET\|CharSet' > $TMPFLD/unsplit_charsets
    # Unassigned1
    # Unassigned2
    # ERROR when regions overlap
}

split_matrix_into_partitions()
#   This function takes the matrix of a NEXUS file and extracts from it a 
#   character range (hereafter "sub-matrix") specified by the definition 
#   of a charsets line.
#   INP:  $1: path to temp folder
#         $2: complete NEXUS file
#         $3: unsplit_matrix
#         $4: unsplit_charsets
#   OUP:  file with name $CHARRNG_FNAME in temp folder
#         file with name $PARTITION_FNAME in temp folder
{
    NEXUS_P1=$(sed -e '/matrix\|MATRIX\|Matrix/,$d' $2)
    NEXUS_P2='\nMATRIX'
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
        echo -e "$NEXUS_P1 $NEXUS_P2" > $1/$PARTITION_FNAME             # Step 1 of assembling the new partition file
        awk 'NR==FNR{start=$1;lgth=$2-$1+1;next} {print $1, substr($2,start,lgth)}' FS='-' $1/$CHARRNG_FNAME FS=' ' $1/$3 >> $1/$PARTITION_FNAME  # Step 2 of assembling the new partition file: Add the sub-matrix
        echo -e "$NEXUS_P4" >> $1/$PARTITION_FNAME                      # Step 3 of assembling the new partition file
        NEW_LENGTH=$(awk '{print $2-$1+1}' FS='-' $TMPFLD/$CHARRNG_FNAME)   # Get the length of the sub-matrix
        sed -i "/dimensions\|DIMENSIONS\|Dimensions/ s/NCHAR\=.*\;/NCHAR\=$NEW_LENGTH\;/" $1/$PARTITION_FNAME  # Replace the number of characters with the length of the sub-matrix
    done < $1/$4
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
echo -ne " INFO  | $(CURRENTTIME) | Step 01: Checking infiles\n"

# Checking if input files exists
if [[ ! -f $NEXUS ]]; then 
    echo -ne " ERROR | $(CURRENTTIME) | File not found: $NEXUS\n"
    exit 1
fi
if [[ ! -f $MODELTEST ]]; then 
    echo -ne " ERROR | $(CURRENTTIME) | File not found: $MODELTEST\n"
    exit 1
fi

## Checking if input files are specified as absolute paths
#if [[ ! $NEXUS = /* ]]; then
#    echo -ne " ERROR | $(CURRENTTIME) | File path not absolute: $NEXUS\n"
#    exit 1
#fi
#if [[ ! $MODELTEST = /* ]]; then
#    echo -ne " ERROR | $(CURRENTTIME) | File path not absolute: $MODELTEST\n"
#    exit 1
#fi

# Define outfile namestem
BASENAME=$(basename $NEXUS) # Using basename to strip off path
FILESTEM=${BASENAME%.nex*} # Using parameter expansion to remove file extension

# Define outfile names
DATABLOC=${FILESTEM}_DataBlock
SETSBLOC=${FILESTEM}_SetsBlock
MDLOVRVW=${FILESTEM}_ModelOverview
LSETSPEC=${FILESTEM}_LsetSpecs
OUTFNAME=${FILESTEM}_partitioned.mrbayes

# Parsing charsets from NEXUS file
CHARSETS=$(cat $NEXUS | grep 'charset\|CharSet\|CHARSET')

# Evaluate number of cores available on Desktop machine
CORES=$(grep -c ^processor /proc/cpuinfo)

# Make and cd into temporary folder
TMPFLD=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
mkdir -p $TMPFLD

########################################################################

## Step 02: Splitting NEXUS file into blocks and sub-blocks
echo -ne " INFO  | $(CURRENTTIME) | Step 02: Splitting NEXUS file into blocks and sub-blocks\n"
split_nexus_into_blocks $NEXUS

########################################################################

## Step 03: Splitting matrix into partitions
echo -ne " INFO  | $(CURRENTTIME) | Step 03: Splitting matrix into partitions\n"
split_matrix_into_partitions $TMPFLD $NEXUS unsplit_matrix unsplit_charsets

########################################################################

## Step 04: Conducting modeltesting via jModelTest2
echo -ne " INFO  | $(CURRENTTIME) | Step 04: Conducting modeltesting via jModelTest2\n"

count=`ls -1 $TMPFLD/partition_* 2>/dev/null | wc -l`
if [[ $count != 0 ]];
then
    for file in ./$TMPFLD/partition_*; do 
    java -jar $MODELTEST -tr $CORES -d "$file" -g 4 -i -f -AIC -o ${file}.bestModel 1>${file}.bestModel.log 2>&1
    done
else
    echo -ne " ERROR | $(CURRENTTIME) | No partition files (partition_*) found.\n"
    exit 1
fi

########################################################################

## Step 05: Extracting model information from jModelTest2
echo -ne " INFO  | $(CURRENTTIME) | Step 05: Extracting model information from jModelTest2\n"

count=`ls -1 $TMPFLD/*.bestModel 2>/dev/null | wc -l`
if [[ $count != 0 ]];
then
    for file in ./$TMPFLD/*.bestModel; do 
    echo -ne "$file" | sed 's/partition_//g' | sed 's/\.bestModel//g' >> $TMPFLD/$MDLOVRVW
    cat "$file" | grep -A1 ' Model selected:' | tail -n1 | sed 's/Model = //g' | sed 's/   / /g' >> $TMPFLD/$MDLOVRVW  # Have to use command ´cat´, cannot use ´echo´ here; not sure why.
    done
else
    echo -ne " ERROR | $(CURRENTTIME) | No result files of the modeltesting (*.bestModel) found.\n"
    exit 1
fi

########################################################################

## Step 06: Converting modeltest results to lset specs for MrBayes
echo -ne " INFO  | $(CURRENTTIME) | Step 06: Converting modeltest results to lset specs for MrBayes\n"

if [[ ! -f $TMPFLD/$MDLOVRVW ]];
then 
    echo -ne " ERROR | $(CURRENTTIME) | File not found: $TMPFLD/$MDLOVRVW\n"
    exit 1
fi

awk -F'/' '{print "[# "$NF}' $TMPFLD/$MDLOVRVW > $TMPFLD/$LSETSPEC

# Retrieving the list of best fitting models at this point
BEST_MDLS=$(sed 's/\[\# //g' $TMPFLD/$LSETSPEC)

# Setting up lset specifications
convert_models_into_lset $TMPFLD $LSETSPEC

########################################################################

## Step 07: Replacing partition names with partition numbers (This is why ordered partitions are CRITICAL!)
echo -ne " INFO  | $(CURRENTTIME) | Step 07: Replacing partition names with partition numbers\n"

# Replacing the partition name sin $LSETSPEC with the line numbers in the NEXUS file
for file in $(echo "$CHARSETS" | awk '{print $2}'); do # Use doublequote-enclosing around a variable to maintain its newlines!
LINENUM=$(echo "$CHARSETS" | awk '/'$file'/{print NR; exit}');
sed -i "/$file/ s/applyto\=()/applyto\=($LINENUM)/" $TMPFLD/$LSETSPEC # Double quotes are critical here due to variable
done

########################################################################

## Step 08: Assembling output file
echo -ne " INFO  | $(CURRENTTIME) | Step 08: Assembling output file\n"

if [[ $OUTFNAME = /* ]];
then
    echo -ne " ERROR | $(CURRENTTIME) | Outfile already exists in filepath: $OUTFNAME\n"
    exit 1
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
echo -ne " INFO  | $(CURRENTTIME) | Step 09: Generating MrBayes block\n"

write_mrbayes_block $TMPFLD $OUTFNAME $LSETSPEC "$CHARSETS"

########################################################################

## Step 10: Cleaning up the work directory
echo -ne " INFO  | $(CURRENTTIME) | Step 10: Cleaning up\n"

rm -r $TMPFLD

echo -e " INFO  | $(CURRENTTIME) | Processing complete for: $NEXUS\n" # End of file message
exit 0 # exit without mistakes

########################################################################
