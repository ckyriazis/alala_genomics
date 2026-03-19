

cd /home/centos/USS/chris/alala/output/data_processing/03_markDups

for bam in *.bam; do
  sample=$(basename "$bam" _markDup.sorted.bam)
  mosdepth \
    --by 100 \
    --chrom NC_063222.1 \
    "$sample" \
    "$bam"
done
