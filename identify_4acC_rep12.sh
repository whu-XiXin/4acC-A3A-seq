#!/bin/bash


#bismark合并文件过滤输出
awk -F'\t' -v low_rate=0.05 -v high_rate=0.15 '
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i;
    print;
    next;
}
$(col["CtoT-rep1_N"])+0 >= 100 &&
$(col["CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_N"])+0 >= 100 &&
$(col["ac4CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_Methylation_rate"])+0 < low_rate &&
$(col["ac4CtoT-rep2_Methylation_rate"])+0 < low_rate &&
$(col["CtoT-rep1_Methylation_rate"])+0 > high_rate 
' /home/bismark/4acC_A3A_seq_combine.txt > /home/bismark/4acC_A3A_seq_combine.Csite.rep1.txt

awk -F'\t' -v low_rate=0.05 -v high_rate=0.15 '
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i;
    print;
    next;
}
$(col["CtoT-rep1_N"])+0 >= 100 &&
$(col["CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_N"])+0 >= 100 &&
$(col["ac4CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_Methylation_rate"])+0 < low_rate &&
$(col["ac4CtoT-rep2_Methylation_rate"])+0 < low_rate &&
$(col["CtoT-rep2_Methylation_rate"])+0 > high_rate
' /home/bismark/4acC_A3A_seq_combine.txt > /home/bismark/4acC_A3A_seq_combine.Csite.rep2.txt

awk -F'\t' -v low_rate=0.05 -v high_rate=0.15 '
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i;
    print;
    next;
}
$(col["CtoT-rep1_N"])+0 >= 100 &&
$(col["CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_N"])+0 >= 100 &&
$(col["ac4CtoT-rep2_N"])+0 >= 100 &&
$(col["ac4CtoT-rep1_Methylation_rate"])+0 < low_rate &&
$(col["ac4CtoT-rep2_Methylation_rate"])+0 < low_rate &&
$(col["CtoT-rep2_Methylation_rate"])+0 > high_rate
$(col["CtoT-rep1_Methylation_rate"])+0 > high_rate
' /home/bismark/4acC_A3A_seq_combine.txt > /home/bismark/4acC_A3A_seq_combine.Csite.txt

