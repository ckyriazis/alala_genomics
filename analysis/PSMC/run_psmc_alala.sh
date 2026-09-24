# PSMC quick guide

# Requirements:
#	- VCF file that has been filtered (genotypes failing filters set to missing)
#	- A tab-delimited file with <chromosome name> in 1st column and <length> in 2nd column ($GENOME below)
#	- Reference genome fasta
#	- Software: bcftools, bedtools, bgzip (can be replaced with gzip?), psmc, gnuplot

# Recommended: simplify VCF files beforehand to reduce file size and exclude unwanted regions/sites
#	- Clear INFO column and unnecessary FORMAT fields (bcftools annotate -x INFO,FORMAT)
#	- Exclude sites that don't pass filters (bcftools view -f .,PASS)
#	- Concatenate chromosomes/scaffold VCF files in same order as reference genome fasta, excluding small scaffolds and sex chromosomes (bcftools concat)


source /home/centos/USS/conservation-genetics/PPM/anaconda3/etc/profile.d/conda.sh
conda activate annot

### Set variables

GENOME=alala_chroms.txt
VCF=175AL_FilterB_Round1_autosomes_PASS.vcf.gz
VCFdir=/home/centos/USS/chris/alala/output/data_processing/08_mergeVCF
SAMPLE=AL99
REF=/home/centos/USS/chris/alala/reference_genome/GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna


# Mutation rate (per site per generation) and generation time (years) for plotting
MU=5e-9
GEN=6


### Generate a mask consisting of all uncalled sites (get called sites first, then get the complement)
# Time: ~2 hr (nearly all time taken by first step)

BED_CALLED=${VCF%.vcf.gz}_${SAMPLE}_sites_called.bed
BED_UNCALLED=${VCF%.vcf.gz}_${SAMPLE}_sites_uncalled.bed

bcftools view --no-header --exclude-uncalled -s ${SAMPLE} -Ov ${VCFdir}/${VCF} \
| awk '{printf "%s\t%s\t%s\n", $1, $2-1, $2}' \
| bedtools merge -i stdin \
| bgzip > ${BED_CALLED}.gz

zcat ${BED_CALLED}.gz \
| bedtools sort -i stdin -g ${GENOME} \
| bedtools complement -i stdin -g ${GENOME} > ${BED_UNCALLED}


### Generate consensus sequence
# Time: ~1 hr

FASTA=${VCF%.vcf.gz}_${SAMPLE}.fa.gz

bcftools consensus -f ${REF} --iupac-codes --missing N --mask ${BED_UNCALLED} -s ${SAMPLE} ${VCFdir}/${VCF} \
| gzip > ${FASTA}


### Generate PSMC input file
# Time: <5 min

PSMCDIR=/home/centos/USS/chris/alala/scripts/analyses/PSMC/psmc-master
PSMCFA=${FASTA%.fa.gz}.psmcfa
${PSMCDIR}/utils/fq2psmcfa ${FASTA} > ${PSMCFA}


### Run PSMC
# Most commonly used time pattern is -p 4+25*2+4+6
# Time: ~2 hr

PSMC_OUT=${PSMCFA%.psmcfa}.psmc

psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o ${PSMC_OUT} ${PSMCFA}


### Check that likelihood scores rapidly plateau
# Time: <1 min

grep "^LK" ${PSMC_OUT} \
| tail -n +2 \
| cut -f2 \
| gnuplot -p -e "set terminal dumb size 120, 30; set autoscale; plot '-' using 1 with lines notitle"


### Check to see if any intervals after Round 20 have <10 recombinations
# Time: <1 min

tail -n +$(grep -P -n "^RD\t21" ${PSMC_OUT} \
| cut -d':' -f1) ${PSMC_OUT} \
| grep "^RS" \
| awk '$5<10'


### Plot PSMC results (using -R will retain temp files, which can be used to generate your own plots)
# Time: <1 min

PLOT=${PSMC_OUT}.plot

perl ${PSMCDIR}/utils/psmc_plot.pl -R -p -u ${MU} -g ${GEN} ${PLOT} ${PSMC_OUT}

# To specify x- and y-axis limits:
# XMIN=1e3
# XMAX=1e7
# YMAX=50
# perl ${PSMCDIR}/utils/psmc_plot.pl -R -p -x ${XMIN} -X ${XMAX} -Y ${YMAX} -u ${MU} -g ${GEN} ${PLOT} ${PSMC_OUT}



