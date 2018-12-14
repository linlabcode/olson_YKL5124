#!/usr/bin/bash

#cufflinksFolder=/path/to/cufflinks/output
#gtfFile=/path/to/hg19_gtf_withERCC
#bamFolder=/path/to/bams

cd $cufflinksFolder/

echo 'making cuffquant folders'
mkdir A2_5-124_2_1
mkdir A2_5-124_2_2
mkdir A2_5-124_2_3
mkdir A2_5-124_500_1
mkdir A2_5-124_500_2
mkdir A2_5-124_500_3
mkdir A2_531_1
mkdir A2_531_2
mkdir A2_531_3
mkdir A2_531_5-124_2_1
mkdir A2_531_5-124_2_2
mkdir A2_531_5-124_2_3
mkdir A2_531_5-124_500_1
mkdir A2_531_5-124_500_2
mkdir A2_531_5-124_500_3
mkdir A2_DMSO_1
mkdir A2_DMSO_2
mkdir A2_DMSO_3
mkdir A2_THZ1_1
mkdir A2_THZ1_2
mkdir A2_THZ1_3
mkdir H1_5-124_2_1
mkdir H1_5-124_2_2
mkdir H1_5-124_2_3
mkdir H1_5-124_500_1
mkdir H1_5-124_500_2
mkdir H1_5-124_500_3
mkdir H1_531_1
mkdir H1_531_2
mkdir H1_531_3
mkdir H1_531_5-124_2_1
mkdir H1_531_5-124_2_2
mkdir H1_531_5-124_2_3
mkdir H1_531_5-124_500_1
mkdir H1_531_5-124_500_2
mkdir H1_531_5-124_500_3
mkdir H1_DMSO_1
mkdir H1_DMSO_2
mkdir H1_DMSO_3
mkdir H1_THZ1_1
mkdir H1_THZ1_2
mkdir H1_THZ1_3

echo 'calling cuffquant'
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_2_1/ $gtfFile $bamFolder/09148_20180202_A2_5-124_2_1_CO5212-1of2_S10_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_2_2/ $gtfFile $bamFolder/09164_20180202_A2_5-124_2_2_CO5212-1of2_S11_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_2_3/ $gtfFile $bamFolder/09166_20180208_A2_5-124_2_3_CO5212-12_S1_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_500_1/ $gtfFile $bamFolder/09153_20180202_A2_5-124_500_1_CO5212-1of2_S7_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_500_2/ $gtfFile $bamFolder/09178_20180202_A2_5-124_500_2_CO5212-1of2_S8_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_5-124_500_3/ $gtfFile $bamFolder/09157_20180202_A2_5-124_500_3_CO5212-1of2_S9_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_1/ $gtfFile $bamFolder/09168_20180202_A2_531_1_CO5212-1of2_S13_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_2/ $gtfFile $bamFolder/09160_20180202_A2_531_2_CO5212-1of2_S14_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_3/ $gtfFile $bamFolder/09154_20180202_A2_531_3_CO5212-1of2_S15_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_2_1/ $gtfFile $bamFolder/09146_20180202_A2_531_5-124_2_1_CO5212-1of2_S19_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_2_2/ $gtfFile $bamFolder/09162_20180202_A2_531_5-124_2_2_CO5212-1of2_S20_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_2_3/ $gtfFile $bamFolder/09174_20180202_A2_531_5-124_2_3_CO5212-1of2_S21_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_500_1/ $gtfFile $bamFolder/09156_20180202_A2_531_5-124_500_1_CO5212-1of2_S16_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_500_2/ $gtfFile $bamFolder/09147_20180202_A2_531_5-124_500_2_CO5212-1of2_S17_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_531_5-124_500_3/ $gtfFile $bamFolder/09185_20180202_A2_531_5-124_500_3_CO5212-1of2_S18_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_DMSO_1/ $gtfFile $bamFolder/09177_20180202_A2_DMSO_1_CO5212-1of2_S1_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_DMSO_2/ $gtfFile $bamFolder/09182_20180202_A2_DMSO_2_CO5212-1of2_S2_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_DMSO_3/ $gtfFile $bamFolder/09176_20180202_A2_DMSO_3_CO5212-1of2_S3_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_THZ1_1/ $gtfFile $bamFolder/09161_20180202_A2_THZ1_1_CO5212-1of2_S4_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_THZ1_2/ $gtfFile $bamFolder/09189_20180202_A2_THZ1_2_CO5212-1of2_S5_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/A2_THZ1_3/ $gtfFile $bamFolder/09172_20180202_A2_THZ1_3_CO5212-1of2_S6_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_2_1/ $gtfFile $bamFolder/09158_20180202_H1_5-124_2_1_CO5212-2of2_S10_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_2_2/ $gtfFile $bamFolder/09170_20180202_H1_5-124_2_2_CO5212-2of2_S11_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_2_3/ $gtfFile $bamFolder/09173_20180202_H1_5-124_2_3_CO5212-2of2_S12_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_500_1/ $gtfFile $bamFolder/09167_20180202_H1_5-124_500_1_CO5212-2of2_S7_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_500_2/ $gtfFile $bamFolder/09186_20180202_H1_5-124_500_2_CO5212-2of2_S8_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_5-124_500_3/ $gtfFile $bamFolder/09171_20180202_H1_5-124_500_3_CO5212-2of2_S9_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_1/ $gtfFile $bamFolder/09159_20180202_H1_531_1_CO5212-2of2_S13_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_2/ $gtfFile $bamFolder/09150_20180202_H1_531_2_CO5212-2of2_S14_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_3/ $gtfFile $bamFolder/09184_20180202_H1_531_3_CO5212-2of2_S15_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_2_1/ $gtfFile $bamFolder/09169_20180202_H1_531_5-124_2_1_CO5212-2of2_S19_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_2_2/ $gtfFile $bamFolder/09149_20180202_H1_531_5-124_2_2_CO5212-2of2_S20_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_2_3/ $gtfFile $bamFolder/09155_20180202_H1_531_5-124_2_3_CO5212-2of2_S21_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_500_1/ $gtfFile $bamFolder/09152_20180202_H1_531_5-124_500_1_CO5212-2of2_S16_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_500_2/ $gtfFile $bamFolder/09165_20180202_H1_531_5-124_500_2_CO5212-2of2_S17_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_531_5-124_500_3/ $gtfFile $bamFolder/09163_20180202_H1_531_5-124_500_3_CO5212-2of2_S18_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_DMSO_1/ $gtfFile $bamFolder/09179_20180202_H1_DMSO_1_CO5212-2of2_S1_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_DMSO_2/ $gtfFile $bamFolder/09180_20180202_H1_DMSO_2_CO5212-2of2_S2_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_DMSO_3/ $gtfFile $bamFolder/09151_20180202_H1_DMSO_3_CO5212-2of2_S3_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_THZ1_1/ $gtfFile $bamFolder/09181_20180202_H1_THZ1_1_CO5212-2of2_S4_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_THZ1_2/ $gtfFile $bamFolder/09175_20180202_H1_THZ1_2_CO5212-2of2_S5_R1_001.bam --library-type fr-firststrand
# cuffquant -p 4 -o $cufflinksFolder/H1_THZ1_3/ $gtfFile $bamFolder/09183_20180202_H1_THZ1_3_CO5212-2of2_S6_R1_001.bam --library-type fr-firststrand  

echo 'running cuffnorm command'
#cuffnorm -p 16 -o $cufflinksFolder/CDK7_GRAY_cuffnorm/ -L A2_5-124_2,A2_5-124_500,A2_531,A2_531_5-124_2,A2_531_5-124_500,A2_DMSO,A2_THZ1,H1_5-124_2,H1_5-124_500,H1_531,H1_531_5-124_2,H1_531_5-124_500,H1_DMSO,H1_THZ1 $gtfFile $cufflinksFolder/A2_5-124_2_1/abundances.cxb,$cufflinksFolder/A2_5-124_2_2/abundances.cxb,$cufflinksFolder/A2_5-124_2_3/abundances.cxb $cufflinksFolder/A2_5-124_500_1/abundances.cxb,$cufflinksFolder/A2_5-124_500_2/abundances.cxb,$cufflinksFolder/A2_5-124_500_3/abundances.cxb $cufflinksFolder/A2_531_1/abundances.cxb,$cufflinksFolder/A2_531_2/abundances.cxb,$cufflinksFolder/A2_531_3/abundances.cxb $cufflinksFolder/A2_531_5-124_2_1/abundances.cxb,$cufflinksFolder/A2_531_5-124_2_2/abundances.cxb,$cufflinksFolder/A2_531_5-124_2_3/abundances.cxb $cufflinksFolder/A2_531_5-124_500_1/abundances.cxb,$cufflinksFolder/A2_531_5-124_500_2/abundances.cxb,$cufflinksFolder/A2_531_5-124_500_3/abundances.cxb $cufflinksFolder/A2_DMSO_1/abundances.cxb,$cufflinksFolder/A2_DMSO_2/abundances.cxb,$cufflinksFolder/A2_DMSO_3/abundances.cxb $cufflinksFolder/A2_THZ1_1/abundances.cxb,$cufflinksFolder/A2_THZ1_2/abundances.cxb,$cufflinksFolder/A2_THZ1_3/abundances.cxb $cufflinksFolder/H1_5-124_2_1/abundances.cxb,$cufflinksFolder/H1_5-124_2_2/abundances.cxb,$cufflinksFolder/H1_5-124_2_3/abundances.cxb $cufflinksFolder/H1_5-124_500_1/abundances.cxb,$cufflinksFolder/H1_5-124_500_2/abundances.cxb,$cufflinksFolder/H1_5-124_500_3/abundances.cxb $cufflinksFolder/H1_531_1/abundances.cxb,$cufflinksFolder/H1_531_2/abundances.cxb,$cufflinksFolder/H1_531_3/abundances.cxb $cufflinksFolder/H1_531_5-124_2_1/abundances.cxb,$cufflinksFolder/H1_531_5-124_2_2/abundances.cxb,$cufflinksFolder/H1_531_5-124_2_3/abundances.cxb $cufflinksFolder/H1_531_5-124_500_1/abundances.cxb,$cufflinksFolder/H1_531_5-124_500_2/abundances.cxb,$cufflinksFolder/H1_531_5-124_500_3/abundances.cxb $cufflinksFolder/H1_DMSO_1/abundances.cxb,$cufflinksFolder/H1_DMSO_2/abundances.cxb,$cufflinksFolder/H1_DMSO_3/abundances.cxb $cufflinksFolder/H1_THZ1_1/abundances.cxb,$cufflinksFolder/H1_THZ1_2/abundances.cxb,$cufflinksFolder/H1_THZ1_3/abundances.cxb --library-type fr-firststrand

Rscript normalizeRNASeq.R $cufflinksFolder/CDK7_GRAY_cuffnorm/genes.fpkm_table $cufflinksFolder/CDK7_GRAY_cuffnorm/output/ CDK7_GRAY A2_5-124_2,A2_5-124_500,A2_531,A2_531_5-124_2,A2_531_5-124_500,A2_DMSO,A2_THZ1,H1_5-124_2,H1_5-124_500,H1_531,H1_531_5-124_2,H1_531_5-124_500,H1_DMSO,H1_THZ1 TRUE
