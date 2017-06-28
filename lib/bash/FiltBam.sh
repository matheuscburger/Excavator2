#!/bin/bash

################################################
#
# .bam file filtering 
#
################################################

Bam_File=$1
MAPQ=$2
Output_Folder=$3
Program_Folder=$4
Sample_Name=$5
Target_Name=$6
Assembly=$7




if [ ! -d "$Output_Folder/.tmp" ]; then
	mkdir $Output_Folder/.tmp
fi

mychr=$(cat $Program_Folder/data/targets/$Assembly/$Target_Name/*_chromosome.txt)



for i in $mychr
do
		samtools \view $Bam_File "$i" | cut -f2,3,4,5 | awk '$1!=4' | cut -f2,3,4 | awk '$3==0 || $3>="/$MAPQ/"' | perl -lane 'print "$F[0]\t$F[1]"' > $Output_Folder/.tmp/.FilteredBamtmp.000
		R --slave --args $Sample_Name,$Output_Folder,$i,$Program_Folder,$Target_Name,$Assembly < $Program_Folder/lib/R/MakeReadCount.R
done
