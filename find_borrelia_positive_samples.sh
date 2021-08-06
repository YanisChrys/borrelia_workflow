#!/bin/bash

touch false_positive_borrelia.txt

while IFS= read -r sample_pos
do
	for my_folder in data/input/multi_refs/*
	do 
	run_id=$(basename $my_folder) 
		if [ -a ${my_folder}/${sample_pos}.bam ]
		then
		cp data/input/multi_refs/${run_id}/${sample_pos}.bam data/input/multi_refs_positive/${run_id}/
		test `samtools view -c -F 4 data/input/multi_refs_positive/${run_id}/${sample_pos}.bam` -le 500 && echo ${sample_pos} >> false_positive_borrelia.txt
		fi 
	done
done < data/metadata/borrelia_pcr.txt