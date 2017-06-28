#!/bin/bash

allchr=$(cat $3/data/targets/$4/$5/*_chromosome.txt)


############################################################################  	if output files already exist print warning

if [[ -e "$3/data/targets/$4/$5/MAP/*RData" ]]; then 
  echo 'Cleaning old files from  Mappability folder!'
  rm -f $3/data/targets/$4/$5/MAP/*RData
fi

if [[ -e "$3/data/targets/$4/$5/GCC/*RData" ]]; then 
	echo 'Cleaning old files from GC-Content folder!'
	rm -f $3/data/targets/$4/$5/GCC/*RData
fi

############################################################################  	fasta file indexing if needed

if [[ ! -e "$6.fai" ]]; then 
	samtools faidx $6
fi


check=$(echo $allchr | grep "chr" | wc -c)
if (($check == 0)); then

for i in $allchr
do
  
  grep -w $i $2 | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print "chr"$1,$2-1,$3,$4}' > $3/data/targets/$4/$5/.temp.bed
  $3/lib/OtherLibrary/./bigWigAverageOverBed $1 $3/data/targets/$4/$5/.temp.bed $3/data/targets/$4/$5/MAP/Mapout.txt
	R --slave --args $i,$3/data/targets/$4/$5/MAP < $3/lib/R/SaveMap.R
  rm -f $3/data/targets/$4/$5/.temp.bed
  rm -f $3/data/targets/$4/$5/MAP/Mapout.txt
  grep -w $i $2 | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2-1,$3,$4}' > $3/data/targets/$4/$5/.temp.bed
  bedtools nuc -fi $6 -bed $3/data/targets/$4/$5/.temp.bed | cut -f6 | sed '1d' > $3/data/targets/$4/$5/GCC/GCC.txt
  R --slave --args $i,$3/data/targets/$4/$5/GCC < $3/lib/R/SaveGCC.R
  rm -f $3/data/targets/$4/$5/.temp.bed
  rm -f $3/data/targets/$4/$5/GCC/GCC.txt
  grep -w $i $2 | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2,$2+1,$4}' > $3/data/targets/$4/$5/.temp.bed
  fastaFromBed -fi $6 -bed $3/data/targets/$4/$5/.temp.bed -fo $3/data/targets/$4/$5/FRB/FRB.txt
  R --slave --args $i,$3/data/targets/$4/$5/FRB < $3/lib/R/SaveFRB.R
  rm -f $3/data/targets/$4/$5/.temp.bed
  rm -f $3/data/targets/$4/$5/FRB/FRB.txt
done

else

for i in $allchr
do
  
  grep -w $i $2 | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2-1,$3,$4}' > $3/data/targets/$4/$5/.temp.bed
  $3/lib/OtherLibrary/./bigWigAverageOverBed $1 $3/data/targets/$4/$5/.temp.bed $3/data/targets/$4/$5/MAP/Mapout.txt
  R --slave --args $i,$3/data/targets/$4/$5/MAP < $3/lib/R/SaveMap.R
  bedtools nuc -fi $6 -bed $3/data/targets/$4/$5/.temp.bed | cut -f6 | sed '1d' > $3/data/targets/$4/$5/GCC/GCC.txt
  R --slave --args $i,$3/data/targets/$4/$5/GCC < $3/lib/R/SaveGCC.R
  rm -f $3/data/targets/$4/$5/.temp.bed
  rm -f $3/data/targets/$4/$5/GCC/GCC.txt
  rm -f $3/data/targets/$4/$5/MAP/Mapout.txt
  grep -w $i $2 | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2,$2+1,$4}' > $3/data/targets/$4/$5/.temp.bed
  fastaFromBed -fi $6 -bed $3/data/targets/$4/$5/.temp.bed -fo $3/data/targets/$4/$5/FRB/FRB.txt
  R --slave --args $i,$3/data/targets/$4/$5/FRB < $3/lib/R/SaveFRB.R
  rm -f $3/data/targets/$4/$5/.temp.bed
  rm -f $3/data/targets/$4/$5/FRB/FRB.txt
done

fi

