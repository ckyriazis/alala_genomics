#!/bin/bash


REFERENCE=/home/centos/USS/chris/alala/reference_genome/GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna
GATK=~/USS/conservation-genetics/nwr/bin/GenomeAnalysisTK.jar
OUTDIR=/home/centos/USS/chris/alala/output/data_processing/06_TrimAlternates_VariantAnnotator

VCF=175AL_noFilter_TrimAlternates_chrom5


cd ${OUTDIR}

java -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${VCF}.vcf.gz -selectType SNP -o ${OUTDIR}/${VCF}_snps.vcf


java -jar ${GATK} -R ${REFERENCE} -T VariantsToTable -V ${OUTDIR}/${VCF}${CHR}_snps.vcf -F AN -F BaseQRankSum -F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR --out ${OUTDIR}/${VCF}_SNPs_table.txt
