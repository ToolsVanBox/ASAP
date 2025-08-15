#!/bin/bash

MYFOLDER=$1
SAMPLE=$2
SAMPLE_TYPE='tumor'

R1_FASTQS=($(gcloud storage ls ${MYFOLDER}/** | grep -P "${SAMPLE}.*_1.f.*q.gz"))
R1_FASTQS+=($(gcloud storage ls ${MYFOLDER}/** | grep -P "${SAMPLE}.*_R1_.*.f.*q.gz"))
R1_FASTQS+=($(gcloud storage ls ${MYFOLDER}/** | grep -P "${SAMPLE}.*_R1.f.*q.gz"))

#R1_FASTQS=($(find ${MYFOLDER} -iname "*_1.f*q.gz"))
#R1_FASTQS+=($(find ${MYFOLDER} -iname "*_R1_*.f*q.gz"))

echo "sample,fastq_1,fastq_2,bam,bai,sample_type"

for FASTQ_1 in ${R1_FASTQS[@]}; do 
	BAM=""
	BAI=""
	SAMPLE=$( echo $(basename ${FASTQ_1}) | cut -f1 -d'_' | cut -f1 -d'.')
        if [[ "${FASTQ_1}" == *"_R1_"* ]]; then
		FASTQ_2=${FASTQ_1/_R1_/_R2_}
	fi
	if [[ "${FASTQ_1}" == *"_1.f"* ]]; then
		FASTQ_2=${FASTQ_1/_1.f/_2.f}
        fi
        if [[ "${FASTQ_1}" == *"_R1.f"* ]]; then
    		FASTQ_2=${FASTQ_1/_R1.f/_R2.f}
    	fi
	if gsutil -q stat ${FASTQ_2}; then
		echo ${SAMPLE},${FASTQ_1},${FASTQ_2},${BAM},${BAI},${SAMPLE_TYPE}
	else
		echo "NO PAIRED FASTQ FILES FOUND FOR ${SAMPLE}"
        fi
done

BAMS=($(gcloud storage ls ${MYFOLDER}/** | grep -P ".*.bam$"))
#BAMS=($(find ${MYFOLDER} -iname "*bam"))

for BAM in ${BAMS[@]}; do 
	FASTQ_1=""
	FASTQ_2=""
	SAMPLE=$( echo $(basename ${BAM}) | cut -f1 -d'_' | cut -f1 -d'.')
	BAI=${BAM}.bai
	
	if ! gsutil -q stat ${BAI} ; then
		BAI=${BAM/.bam/.bai}
	fi
	if ! gsutil -q stat ${BAI}; then
		echo "NO BAI FILE FOUND FOR ${SAMPLE}"
	else
		echo ${SAMPLE},${FASTQ_1},${FASTQ_2},${BAM},${BAI},${SAMPLE_TYPE}
	fi
done

