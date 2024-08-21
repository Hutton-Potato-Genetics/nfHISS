#!/bin/bash
Reads=$(cat $1 | grep -m 1 'reads' | cut -f5 -d ' ')
Bases=$(cat $1 | grep -m 1 'bases' | cut -f5 -d ' ')
printf "Reads\tBases\n" > $2
printf "$Reads\t$Bases" >> $2
