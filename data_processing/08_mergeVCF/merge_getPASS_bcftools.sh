
indir=/home/centos/USS/chris/alala/output/data_processing/07_FilterVCFfile/
outdir=/home/centos/USS/chris/alala/output/data_processing/08_mergeVCF
prefix='175AL_FilterB_Round1'
suffix='.vcf.gz'


bcftools concat -Oz \
${indir}/${prefix}_chrom1${suffix} \
${indir}/${prefix}_chrom2${suffix} \
${indir}/${prefix}_chrom3${suffix} \
${indir}/${prefix}_chrom4${suffix} \
${indir}/${prefix}_chrom5${suffix} \
${indir}/${prefix}_chrom6${suffix} \
${indir}/${prefix}_chrom7${suffix} \
${indir}/${prefix}_chrom8${suffix} \
${indir}/${prefix}_chrom9${suffix} \
${indir}/${prefix}_chrom10${suffix} \
${indir}/${prefix}_chrom11${suffix} \
${indir}/${prefix}_chrom12${suffix} \
${indir}/${prefix}_chrom13${suffix} \
${indir}/${prefix}_chrom14${suffix} \
${indir}/${prefix}_chrom15${suffix} \
${indir}/${prefix}_chrom16${suffix} \
${indir}/${prefix}_chrom17${suffix} \
${indir}/${prefix}_chrom18${suffix} \
${indir}/${prefix}_chrom19${suffix} \
${indir}/${prefix}_chrom20${suffix} \
${indir}/${prefix}_chrom21${suffix} \
${indir}/${prefix}_chrom22${suffix} \
${indir}/${prefix}_chrom23${suffix} \
${indir}/${prefix}_chrom24${suffix} \
${indir}/${prefix}_chrom25${suffix} \
${indir}/${prefix}_chrom26${suffix} \
${indir}/${prefix}_chrom27${suffix} \
${indir}/${prefix}_chrom28${suffix} \
${indir}/${prefix}_chrom29${suffix} \
${indir}/${prefix}_chrom30${suffix} \
${indir}/${prefix}_chrom31${suffix} \
${indir}/${prefix}_chrom32${suffix} \
${indir}/${prefix}_chrom33${suffix} \
${indir}/${prefix}_chrom34${suffix} \
${indir}/${prefix}_chrom35${suffix} \
${indir}/${prefix}_chrom36${suffix} \
${indir}/${prefix}_chrom37${suffix} \
${indir}/${prefix}_chrom38${suffix} \
${indir}/${prefix}_chrom39${suffix} \
${indir}/${prefix}_chrom40${suffix} \
${indir}/${prefix}_chrom41${suffix} | \
bcftools view -i 'FILTER="PASS"' -Oz -o ${outdir}/${prefix}_autosomes_PASS.vcf.gz


tabix -p vcf ${outdir}/${prefix}_autosomes_PASS.vcf.gz





