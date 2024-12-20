#!/bin/bash

MYFOLDER=$1
SAMPLE_TYPE='tumor'

R1_FASTQS=($(find ${MYFOLDER} -iname "*_1.f*q.gz"))
R1_FASTQS+=($(find ${MYFOLDER} -iname "*_R1_*.f*q.gz"))

echo "sample,fastq_1,fastq_2,bam,bai,cram,crai,sample_type"

for FASTQ_1 in ${R1_FASTQS[@]}; do 
	BAM=""
	BAI=""
	CRAM=""
	CRAI=""
	SAMPLE=$( echo $(basename ${FASTQ_1}) | cut -f1 -d'_')
        if [[ "${FASTQ_1}" == *"_R1_"* ]]; then
		FASTQ_2=${FASTQ_1/_R1_/_R2_}
		echo ${FASTQ_2}
	fi
	if [[ "${FASTQ_1}" == *"_1.f"* ]]; then
		FASTQ_2=${FASTQ_1/_1.f/_2.f}
        fi
	if [[ -f ${FASTQ_1} && ${FASTQ_2} ]]; then
		echo ${SAMPLE},${FASTQ_1},${FASTQ_2},${BAM},${BAI},${CRAM},${CRAI},${SAMPLE_TYPE}
	else
		echo "NO PAIRED FASTQ FILES FOUND FOR ${SAMPLE}"
        fi
done

BAMS=($(find ${MYFOLDER} -iname "*bam"))

for BAM in ${BAMS[@]}; do 
	FASTQ_1=""
	FASTQ_2=""
	CRAM=""
	CRAI=""
	SAMPLE=$( echo $(basename ${BAM}) | cut -f1 -d'_')
	BAI=${BAM}.bai
	if [[ ! -f ${BAI} ]]; then
		BAI=${BAM/.bam/.bai}
	fi
	if [[ ! -f ${BAI} ]]; then
		echo "NO BAI FILE FOUND FOR ${SAMPLE}"
	fi
	echo ${SAMPLE},${FASTQ_1},${FASTQ_2},${BAM},${BAI},${CRAM},${CRAI},${SAMPLE_TYPE}
done

CRAMS=($(find ${MYFOLDER} -iname "*cram"))

for CRAM in ${CRAMS[@]}; do 
	FASTQ_1=""
	FASTQ_2=""
	BAM=""
	BAI=""
	SAMPLE=$( echo $(basename ${CRAM}) | cut -f1 -d'_')
	CRAI=${CRAM}.crai
	if [[ ! -f ${CRAI} ]]; then
		CRAI=${CRAM/.cram/.crai}
	fi
	if [[ ! -f ${CRAI} ]]; then
		echo "NO CRAI FILE FOUND FOR ${SAMPLE}"
	fi
	echo ${SAMPLE},${FASTQ_1},${FASTQ_2},${BAM},${BAI},${CRAM},${CRAI},${SAMPLE_TYPE}
done

