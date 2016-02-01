#!/bin/bash
#
#

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters. This script needs the output directory of Brownie"
    echo "Correct usage: makePlot.sh Brownie's_output_directory"
    echo "The output file is saved in the same directory as the input file in pdf format"
    exit 0
fi



outputDir=$1
covArr=($outputDir'/Cov/*.dat')

for k in ${covArr[@]}
do 
a=$k
b=$a
b=$b'.pdf'
plotName=covPlot.dem
gnuplot -e "inputFilename='$a'" -e "outputFilename='$b'" $plotName 
done
outputDir=$1

Arr=($outputDir'/Statistic/*.dat')
for k in ${Arr[@]}
do 
a=$k
b=$a
b=$b'.pdf'
plotName=componentSizePlot.dem
gnuplot -e "inputFilename='$a'" -e "outputFilename='$b'" $plotName 
done



Arr=($outputDir'/N50/N50_0*.dat')
for k in ${Arr[@]}
do 
a=$k
b=$a
b=$b'.pdf'
plotName=n50Plot.dem
gnuplot -e "inputFilename='$a'" -e "outputFilename='$b'" $plotName 
done





