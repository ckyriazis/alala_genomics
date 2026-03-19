


indir=/home/centos/USS/chris/alala/output/data_processing/07_FilterVCFfile/
#infile=175AL_FilterB_Round1_autosomes_PASS.vcf.gz
infile=175AL_FilterB_Round1_chrom43.vcf.gz
outdir=/home/centos/USS/chris/alala/output/analyses/ROH/
#outfile=175AL_FilterB_Round1_autosomes_PASS.out
outfile=175AL_FilterB_Round1_Zchrom.out

# need -G flag to indicate genotypes rather than genotype likelihoods

bcftools roh -G30 ${indir}${infile} --include 'FILTER="PASS"' -o ${outdir}${outfile}

## this ran for like 2 weeks and did not finish
#bcftools roh -e - -G 30 ${indir}${infile} --include 'FILTER="PASS"' -o ${outdir}${outfile}

