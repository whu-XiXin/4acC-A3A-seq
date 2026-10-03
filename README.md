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

Before running bismark.sh, please install fastp, bismark, bowtie2, samtools and other necessary software on your own. This script will output txt files including C to T rate in different group.

```
bash bismark.sh
```

**Runing combine_data.py**

This script will combine C to T rate in different group and output txt frame.

```
python combine_data.py
```

**Runing identify_4acC_rep12.sh**

This script will identify the 4acC site in C to T rate file.

```
bash identify_4acC_rep12.sh
```

**Runing fisher_test.R**

This script will calculate the fisher test p value of 4acC site

```
Rscript fisher_test.R
```
