#!/bin/bash

SCRPT=/home/michael_science/git/michaelgruenstaeudl_partitionedNex2MrBayes/partitionedNex2MrBayes_2017.09.14.1900.sh
MDLTBIN=/home/michael_science/binaries/jModelTest_2.1.7/jModelTest.jar
IND=/home/michael_science/git/michaelgruenstaeudl_partitionedNex2MrBayes/examples/01_input
INF=test3.nex
FILESTEM=${INF%.nex*}

# WRITE SHELL COMMANDS INTO FILE
#echo '#!/bin/bash' > ${FILESTEM}.shellcmd
#echo "bash $SCRPT -f $IND/$INF -m $MDLTBIN -k -v 1> ${FILESTEM}.log 2>&1" >> ${FILESTEM}.shellcmd
#chmod a+x ${FILESTEM}.shellcmd


# WRAP CMD IN VALGRIND/MASSIF AND EXECUTE CMD
#valgrind --tool=massif --trace-children=yes --pages-as-heap=yes --time-unit=ms --massif-out-file=${FILESTEM}_massif.out.%p ./${FILESTEM}.shellcmd
valgrind --tool=massif --trace-children=yes --pages-as-heap=yes --time-unit=ms --massif-out-file=${FILESTEM}_massif.out.%p bash $SCRPT -f $IND/$INF -m $MDLTBIN -k -v 1> ${FILESTEM}.log 2>&1

# SUMMARIZE VALGRIND/MASSIF OUTPUT
ms_print --threshold=75.0 ${FILESTEM}_massif.out.* > ${FILESTEM}_massif.summary
# For interpretation of graph, see: http://valgrind.org/docs/manual/ms-manual.html

# CLEANING UP
#rm ${FILESTEM}.shellcmd
