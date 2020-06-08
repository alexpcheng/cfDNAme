#!/bin/bash

taxid=$1
taxid_to_gi=$2
ext=$3
curated=$4

awk -v v1=$taxid '$2==v1' $taxid_to_gi | cut -f 1 > temp.gis.$taxid$ext
mkdir -p grefs/$ext/$taxid
blastdbcmd -db $curated -entry_batch temp.gis.$taxid$ext > grefs/$ext/$taxid/genomes.fna
rm temp.gis.$taxid$ext
