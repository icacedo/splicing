#!/bin/bash

FILE1=genome.gz
if test -f "$FILE1"; then
	echo "$FILE1 exists."
else
	ln -s data/c_elegans.PRJNA13758.WS282.genomic.fa.gz genome.gz
fi

FILE2=gff3.gz
if test -f "$FILE2"; then
	echo "$FILE2 exists."
else
	ln -s data/c_elegans.PRJNA13758.WS282.annotations.gff3.gz gff3.gz
fi

FILE3=ws282.gff3
if test -f "$FILE3"; then
	echo "$FILE3 exists."
else
	gunzip -c gff3.gz | grep -E "WormBase|RNASeq" > ws282.gff3
fi

export PYTHONPATH="/home/izzy/Code/grimoire"
export PATH="/home/izzy/Code/grimoire/bin:$PATH"

haman genome.gz ws282.gff3 pcg genes --issuesok

export PYTHONPATH="$PYTHONPATH:/home/izzy/Code/grimoire/grimoire"
