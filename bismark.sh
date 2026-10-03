#!/bin/bash

mkdir ./A3A.Ecoil.rep1

cd ./A3A.Ecoil.rep1

export TMPDIR=./A3A.Ecoil.rep1

bismark -parallel 8 -genome /home/ref/bismark/NCBI_Ecoil/ -1 /home/A3A.Ecoil.rep1/A3A.Ecoil.rep1_1.fq.gz -2 /home/A3A.Ecoil.rep1/A3A.Ecoil.rep1_2.fq.gz --score_min L,0,-0.6 -X 1000 -o /A3A.Ecoil.rep1 --temp_dir /A3A.Ecoil.rep1

deduplicate_bismark --paired --output_dir ./A3A.Ecoil.rep1 ./A3A.Ecoil.rep1/A3A.Ecoil.rep1.bam

samtools sort -@ 9  -l 4 -m 6G -o /A3A.Ecoil.rep1/A3A.Ecoil.rep1.deduplicated.sorted.bam /A3A.Ecoil.rep1/A3A.Ecoil.rep1.deduplicated.bam

mkdir ./A3A.Ecoil.rep1/methylation_data

bismark_methylation_extractor --gzip --cytosine_report --merge_non_CpG --CX_context --comprehensive --genome_folder /home/ref/bismark/NCBI_Ecoil/ \
	./A3A.Ecoil.rep1/A3A.Ecoil.rep1.deduplicated.bam \
	--output_dir ./A3A.Ecoil.rep1/methylation_data --parallel 16

zcat ./A3A.Ecoil.rep1/methylation_data/A3A.Ecoil.rep1.deduplicated.CX_report.txt.gz | awk 'BEGIN {OFS="	"; print "chr", "pos", "strand", "N", "X", "Type", "3-base"} {sum=$4+$5; if (sum >= 5) print $1, $2, $3, sum, $4, $6, $7}'  > ./A3A.Ecoil.rep1/methylation_data/A3A.Ecoil.rep1.deduplicated.CX_report.dss.txt