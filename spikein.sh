#!/bin/bash

mkdir /home/hanjingyu/project/data/postdata/pileup/H8_L3_BS_HeLa_20251018_zf

cd ./

bowtie2 -p 12 -x /home/ref/gDNA_spikein_3base/gDNA_spikein_zf_3base/gDNA_spikein_zf_3base -1 /home/fastp/H8/H8_1.fq.gz -2 /home/fastp/H8/H8_2.fq.gz \
	 -S /home/bowtie/H8.sam 2> /home/bowtie/H8_result.log

samtools view -@ 8 -b /home/bowtie/H8/H8.sam > /home/bowtie/H8/H8.bam

samtools sort -@ 4  -l 9 -m 4G -o /home/bowtie/H8/H8.sorted.bam /home/bowtie/H8/H8.bam

samtools index /home/bowtie/H8/H8.sorted.bam

/home/hanjingyu/project/miniforge3/bin/samtools mpileup -d 8000 -r Spikein \
	--reference /home/ref/gDNA_spikein_3base/gDNA_spikein_zf_3base/gDNA_spikein_zf_3base/gDNA_spikein_zf_3base.fa \
	/home/bowtie/H8/H8.sorted.bam  > /home/bowtie/H8/H8.spikein.pileup

varscan pileup2cns  /home/bowtie/H8/H8.spikein.pileup \
	--min-coverage 5 --min-reads2 1 --min-avg-qual 0 --min-var-freq 15 \
	>  /home/bowtie/H8/H8.spikein.pileup.vcf

/home/hanjingyu/project/miniforge3/bin/samtools mpileup -d 8000 -r Spikein_R \
	--reference /home/ref/gDNA_spikein_3base/gDNA_spikein_zf_3base/gDNA_spikein_zf_3base/gDNA_spikein_zf_3base.fa \
	/home/bowtie/H8/H8.sorted.bam  > /home/bowtie/H8/H8.spikein_R.pileup

varscan pileup2cns  /home/bowtie/H8/H8.spikein_R.pileup \
	--min-coverage 5 --min-reads2 1 --min-avg-qual 0 --min-var-freq 15 \
	>  /home/bowtie/H8/H8.spikein_R.pileup.vcf
