#!/bin/bash

sed -n -e '/KPYR HUMAN/{p;n;}' bipartion.dat > human_bipart.dat
awk -F\' '{ for (i=3; i<=NF; i+=2) $i=""} 1' human_bipart.dat > human_newicks1.dat
cut -d " " -f 4- human_newicks1.dat > human_newicks.dat
awk '{printf "%d\t%s\n", NR, $0}' human_newicks.dat >  human_newicks.dat

mkdir bipartitions
while read -r filename content ; do
	    printf '%s\n' "$content" >> bipartitions/"${filename}.seq"
    done < human_newicks.dat

sed 's/>//' P1-initial.fasta > alignment_file.fasta
