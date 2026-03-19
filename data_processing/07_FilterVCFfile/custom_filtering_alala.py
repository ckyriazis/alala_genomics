'''
Custom filtering of VCF files

Input: raw VCF
Output: filtered VCF prints to screen
- Sites failing filters are marked as FAIL_? or WARN_? in the 7th column
- Sites where REF is in [A,C,G,T] and ALT is in [A,C,G,T,.] go on to genotype filtering if AD and DP present in FORMAT
- Filtered out genotypes are changed to './.', all others reported
- Sites with non-reference alleles remaining after genotype filtering also filtered based on values in INFO column

Possible usage:

SCRIPT=filterVCF.py
python ${SCRIPT} myfile.vcf.gz | bgzip > myfile_filtered.vcf.gz
tabix -p vcf myfile_filtered.vcf.gz

'''

import sys
import gzip
import re

vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'rt')


# Min depth
minD=3

# Max depth (2.5x mean from scaff1)
maxD={
'AL101':45.7,
'AL102':54.13,
'AL103':53.08,
'AL114':60.35,
'AL119':37.4,
'AL120':44.35,
'AL121':38.98,
'AL124':42.78,
'AL125':62.6,
'AL127':64.53,
'AL129':61.35,
'AL133':50.25,
'AL134':18.58,
'AL137':12.1,
'AL138':59.55,
'AL141':64.78,
'AL142':34.63,
'AL143':42.18,
'AL146':62.85,
'AL149':22.2,
'AL153':17.55,
'AL155':67.15,
'AL156':14.63,
'AL159':42,
'AL160':72.9,
'AL162':43.95,
'AL166':37.48,
'AL170':49,
'AL173':60.63,
'AL174':53.88,
'AL175':39.68,
'AL176':64.85,
'AL177':48.5,
'AL178':18.98,
'AL180':41.6,
'AL181':20.75,
'AL182':40.3,
'AL184':42.53,
'AL185':44.38,
'AL187':56.75,
'AL189':38.25,
'AL193':39.28,
'AL196':37.43,
'AL197':46.35,
'AL199':38.15,
'AL200':41.83,
'AL201':38.15,
'AL204':60.23,
'AL207':36.1,
'AL208':13.83,
'AL209':39.8,
'AL210':60.28,
'AL213':42.2,
'AL215':53.68,
'AL216':40.18,
'AL220':60.28,
'AL223':36.58,
'AL231':44.93,
'AL238':55.63,
'AL241':18.8,
'AL242':13.3,
'AL243':45.6,
'AL244':36.08,
'AL246':40.33,
'AL254':53.95,
'AL255':14.65,
'AL256':17.05,
'AL263':74.85,
'AL266':46.45,
'AL269':47.5,
'AL272':57.98,
'AL273':16.25,
'AL274':46.3,
'AL276':45.9,
'AL278':14.48,
'AL279':96.6,
'AL283':72.65,
'AL284':71.3,
'AL288':56.95,
'AL292':46.73,
'AL300':38.28,
'AL303':26.6,
'AL304':44.25,
'AL307':43.43,
'AL309':42.58,
'AL313':63.88,
'AL32':51.38,
'AL35':36.13,
'AL75':50.1,
'AL82':50.08,
'AL86':60.73,
'AL88':50.13,
'AL91':53.23,
'AL92':63.23,
'AL95':54.9,
'AL97':50.45,
'AL99':51.83,
'K1423':40.75,
'K1461':54.28,
'K1542':46.43,
'K1544':37.58,
'K1567':59.98,
'K1579':115.65,
'K1580':41.7,
'K1583':47.88,
'K1603':52.13,
'K1606':68.83,
'K1608':70.33,
'K16100':57.03,
'K1614':50.35,
'K1617':67.25,
'K1618':61.13,
'K1623':20.88,
'K1627':28.28,
'K1632':51.08,
'K1640':75.28,
'K1650':50.58,
'K1660':41.25,
'K1663':40.5,
'K1665':62.43,
'K1669':78.28,
'K1672':81.33,
'K1680':47,
'K1684':76,
'K1693':38.78,
'K1696':44.38,
'K17007':43.43,
'K17008':19.5,
'K17010':66.95,
'K17013':39.45,
'K17014':44.83,
'K17025':52.73,
'K17028':43.08,
'K17030':45.28,
'K17032':110.75,
'K17033':85.28,
'K17034':40.6,
'K17035':34.98,
'K17037':42.33,
'K17040':69.53,
'K17042':60.08,
'K17044':65.23,
'K17045':237.2,
'K17048':41.55,
'K17052':89.03,
'K17053':47.6,
'K17056':75.9,
'K17062':52.9,
'K17064':45.23,
'K17073':45.7,
'K17077':45.38,
'K17078':49.8,
'K17081':54.45,
'K17083':52.5,
'K17087':66.08,
'K17091':50.38,
'K17094':64.5,
'K17098':49.95,
'K17099':48.53,
'K17101':71.95,
'K17102':56.75,
'K17103':66.27,
'K17104':42.88,
'K17106':63.2,
'K18010':45.7,
'K18017':65.28,
'K18019':58.33,
'K18022':48.65,
'K18034':33.63,
'K18078':53.88,
'K18083':49.03,
'K18084':45.75,
'K18091':55.95,
'K19018':50.5,
'K19019':47.88
}

# Individual genotype filtering function
#  - Genotypes failing filters are set to missing (./.)
#  - Applies individual min and max depth filters
#  - Filters heterozygotes if the allele balance (REF/DP) is <20% or >80%
#  - Filters homozygotes if more than 10% of the alleles are different type
#  - 'sample' is the sample name
#  - 'GT_entry' is the entire genotype entry for that individual (typically GT:AD:DP:GQ)
#  - 'ADpos' is the position of the AD field in FORMAT (determined below)
#  - 'DPpos' is the position of the DP field in FORMAT (determined below)

def GTfilter(sample, GT_entry, ADpos, DPpos):
    if GT_entry[:1]=='.' : return GT_entry
    else:
        gt=GT_entry.split(':')
        if gt[0] in ('0/0','0/1','1/1') and gt[DPpos]!='.':
            DP=int(gt[DPpos])
            if minD<=DP<=maxD[sample]:
                REF=float(gt[ADpos].split(',')[0])
                AB=float(REF/DP)
                if gt[0]=='0/0':
                    if AB>=0.9: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                elif gt[0]=='0/1':
                    if 0.2<=AB<=0.8: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                elif gt[0]=='1/1':
                    if AB<=0.1: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                else: './.:' + ':'.join(gt[1:])
            else: return './.:' + ':'.join(gt[1:])
        else: return './.:' + ':'.join(gt[1:])


# Get list of samples in VCF file
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Go back to beginning of file
VCF.seek(0)


# Write pre-existing header lines & add new lines describing filters being applied
for line0 in VCF:
    if line0.startswith('#'):
        if line0.startswith('##FORMAT'):
            sys.stdout.write('##FILTER=<ID=FAIL_REF,Description="Reference allele not one of [A,C,G,T].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ALT,Description="Alternate allele not one of [A,C,G,T,.].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="AD not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="DP not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGT,Description="No called genotypes remain after filtering.">\n')
            sys.stdout.write('##FILTER=<ID=WARN_missing,Description="Excess missingness (>25% of samples uncalled or set to missing).">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_QD,Description="QD < 4.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_FS,Description="FS > 12.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQ,Description="MQ < 40.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQRankSum,Description="MQRankSum < -12.5.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ReadPosRankSum,Description="ReadPosRankSum < -8.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_SOR,Description="SOR > 3.5.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_excessHet,Description="Excess heterozygosity (>75% of genotypes are 0/1).">\n')
            sys.stdout.write(line0)
            break
        else: sys.stdout.write(line0)


# Go through VCF file line by line to apply filters
for line0 in VCF:
    if line0.startswith('#'):
        sys.stdout.write(line0)
        continue
    line=line0.strip().split('\t')

### Site filtering:
### Keep any filters that have already been applied
    filter=[]
    if line[6] not in ('.', 'PASS'):
        filter.append(line[6])

### Check REF allele
    if line[3] not in ['A','C','G','T']:       
        filter.append('FAIL_REF') 

### Check ALT allele
    if line[4] not in ['A','C','G','T','.']:
        filter.append('FAIL_ALT') 

### Access INFO field annotations
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        if '=' in INFO:
            d={INFO.split('=')[0]:INFO.split('=')[1]}
        else: d={}

### Get the position of AD, DP in genotype fields
    if 'AD' in line[8]:
        ADpos=line[8].split(':').index('AD')
    else: filter.append('FAIL_noADi')

    if 'DP' in line[8]:
        DPpos=line[8].split(':').index('DP')
    else: filter.append('FAIL_noDPi')

### If any filters failed, write out line and continue
    if filter!=[]:
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) )
        continue

### Genotype filtering:
    GT_list=[]
    for i in range(0,len(samples)):
        GT=GTfilter(samples[i],line[i+9],ADpos,DPpos)
        GT_list.append(GT)

### Recalculate AC, AN, AF for INFO (after this step, modified INFO values will be output)
    REF=2*[x[:3] for x in GT_list].count('0/0') + [x[:3] for x in GT_list].count('0/1')
    ALT=2*[x[:3] for x in GT_list].count('1/1') + [x[:3] for x in GT_list].count('0/1')
    if REF+ALT==0:
        filter.append('FAIL_noGT')
        sys.stdout.write('%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:9]), '\t'.join(GT_list)) )
        continue    
    d['AC']=ALT
    d['AN']=REF+ALT
    d['AF']=round(float(ALT)/(float(REF)+float(ALT)), 4)

### Warn if >25% of genotypes missing
    n_missing=sum(x[:3]=='./.' for x in GT_list)
    if n_missing>0.25*len(samples):
        filter.append('WARN_missing')

### Fail sites with excess heterozygosity (>75% of genotypes are heterozygous)
    n_het=sum(x[:3]=='0/1' for x in GT_list)
    n_called=sum(x[:3]!='./.' for x in GT_list)
    if float(n_het/n_called)>0.75:
        filter.append('FAIL_excessHet')
        
### Set VariantType, outputting sites with just hom. REF genotypes without further filtering
    if ALT==0:
        d['VariantType']='NO_VARIATION'   
        if filter==[]:
            filter.append('PASS')
        sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )
        continue
    elif REF==0:
        d['VariantType']='NO_VARIATION'
    else:
        d['VariantType']='SNP'

### Fail sites with poor variant metrics
    if 'QD' in d and float(d['QD']) < 4.0:
        filter.append('FAIL_QD')
    if 'FS' in d and float(d['FS']) > 12.0:
        filter.append('FAIL_FS')
    if 'MQ' in d and float(d['MQ']) < 40.0:
        filter.append('FAIL_MQ')
    if 'MQRankSum' in d and float(d['MQRankSum']) < -12.5:
        filter.append('FAIL_MQRankSum')
    if 'ReadPosRankSum' in d and float(d['ReadPosRankSum']) < -8.0:
        filter.append('FAIL_ReadPosRankSum')
    if 'SOR' in d and float(d['SOR']) > 3:
        filter.append('FAIL_SOR')

### Write out new line
    if filter==[]:
        filter.append('PASS')
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )


# Close files and exit
VCF.close()
exit()


