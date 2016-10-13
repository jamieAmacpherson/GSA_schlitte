#! /usr/bin/bash

sed -n -e '/KPYR HUMAN/{p;n;}' bipartion.dat > human_bipart.dat
awk -F\' '{ for (i=3; i<=NF; i+=2) $i=""} 1' human_bipart.dat > human_newicks1.dat
cut -d " " -f 4- human_newicks1.dat > human_newicks.dat
