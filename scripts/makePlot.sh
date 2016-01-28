#!/bin/bash
#
#

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters. This script needs two essential parameters"
    echo "Correct usage: makePlot.sh inputFilename plotName"
    echo "The output file is saved in the same directory as the input file in pdf format"
    exit 0
fi

a=$1 
b=$a
b=$b'.pdf'
plotName=$2
gnuplot -e "inputFilename='$a'" -e "outputFilename='$b'" $plotName 