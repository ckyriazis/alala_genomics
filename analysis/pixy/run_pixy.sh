

indir=/home/centos/USS/chris/alala/output/data_processing/08_mergeVCF
vcf=174AL_FilterB_Round2_autosomes_PASS.vcf.gz

pixy --stats pi \
--vcf ${indir}/${vcf} \
--populations popfile_alive.txt \
--window_size 100000 \
--n_cores 16 \
--output_folder /home/centos/USS/chris/alala/output/analyses/pixy

