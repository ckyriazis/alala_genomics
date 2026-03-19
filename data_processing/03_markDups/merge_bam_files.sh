


REFERENCE=/home/centos/USS/chris/hawaiian_birds/reference_genomes/akikiki-hap2/akikiki-san3754-hap2.fa

indir=/home/centos/USS/chris/hawaiian_birds/output/data_processing/03_markDups



picard MergeSamFiles \
      TMP_DIR=${indir} \
      I=${indir}/AKIK_AI002_redo_markDup.sorted.bam  \
      I=${indir}/AKIK_AI002_markDup.sorted.bam \
      O=${indir}/AKIK_AI002_markDup_merged.sorted.bam
      R=$REFERENCE

picard MergeSamFiles \
      TMP_DIR=${indir} \
      I=${indir}/AKIK_AI043_redo_markDup.sorted.bam  \
      I=${indir}/AKIK_AI043_markDup.sorted.bam \
      O=${indir}/AKIK_AI043_markDup_merged.sorted.bam      
      R=$REFERENCE



