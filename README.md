# 4acC-A3A-seq
----------------------------------------
## Scripts for analysing 4acC-A3A-seq data ##
----------------------------------------
Tools for analyzing data from the 4acC-A3A-seq to identify 4acC modified in spikein and genome.
----------------------------------------
### The link address:
Github: https://github.com/whu-XiXin/4acC-A3A-seq/
-----------------------------------------

## Data analysis process
------------------------------------	

### For spikein in gDNA data 


**Running spikein.sh**

Before running spikein.sh, please install fastp, bowtie2, varscan, samtools and other necessary software on your own. This script will output spike mutation vcf files both forward and reverse strands.

```
bash spikein.sh
```
### For genome data

**Runing bismark.sh**

Before running bismark.sh, please install fastp, cutadapt, bismark, bowtie2, samtools and other necessary software on your own. This script will output txt files including C to T rate in different group.

```
bash bismark.sh
```

**Runing 4acC_identify.R**

This script will calculate the 4acC site

```
Rscript 4acC_identify.R
```
