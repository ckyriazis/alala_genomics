


indir=/home/centos/USS/chris/alala/output/data_processing/08_mergeVCF
infile=175AL_FilterB_Round1_autosomes_PASS
outdir=/home/centos/USS/chris/alala/output/data_processing/09_getPassSNPs


bcftools view \
    --types snps \
    -m2 -M2 \
    --include 'FILTER="PASS" && AC>0' \
    -Oz \
    -o ${outdir}/${infile}_SNPs.vcf.gz \
    ${indir}/${infile}.vcf.gz

tabix -p vcf ${outdir}/${infile}_SNPs.vcf.gz


