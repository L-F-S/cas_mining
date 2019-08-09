#!/bin/bash
# Fri Jul 5 16:31:15 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
if [ $# -lt 1 ]; then
echo ++++++++++++++++++++++++++++++++++++++++++++
echo Usage: 
echo    ./retrieve_minced_annotations_of_dataset\.sh \<dataset\>
exit 1
fi
					
dataset=$1
echo $dataset
cd /shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno
for i in *; do
	if [[ $i == ${dataset}* ]]; then 
		echo $i
		less $i | grep CRISPR | grep minced >/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/justminced/${dataset}/${i}.minced
       	fi; done
