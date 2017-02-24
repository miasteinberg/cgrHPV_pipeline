# vim: ft=python
import os
import glob
import itertools
import pandas
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

configfile: "hpv_config.yaml"
workdir: os.environ['PWD']
shell.executable('bash')

PANEL = config['panel_type']
hpv_type = config['hpv_type']
if PANEL == 'single':
    hpv_ref = config['hpv_reference']['single'] %(hpv_type, hpv_type)
    amp_bed = config['amplicon_bed']['single'] %hpv_type
    len_bed = config['len_bed']['single'] %(hpv_type, hpv_type)
    annot = config['annotation']['single'] %hpv_type
else:
    hpv_ref = config['hpv_reference'][PANEL]
    amp_bed = config['amplicon_bed'][PANEL]
    len_bed = config['len_bed'][PANEL]
    annot = config['annotation'][PANEL]

gatk = config["gatk"]
cov_dev = config["cov_dev"]
GENES = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'L1', 'L2']

bed = {}
len_file = open(len_bed, 'r')
# create a dictionary of type and genome length
# e.g. bed['HPV16_Ref'] = 7906
for line in len_file:
    (type, length) = (line.split()[0], int(line.split()[2]))
    if PANEL == 'hpv16':
        bed[type] = length
    else:
        bed[type] = length - 400
len_file.close()

## ---- The parser will have to be customized for each run ---- ##
def parse_sampleID(filename):
    base = filename.split('/')[-1] 
    # For filenames IonXpress_046_SC075985.tmap.bam
    if base.startswith('I'):
        return filename.split('/')[-1].split('_')[2].split('.')[0]
    # otherwise PAP1111_2001_IonXpress_20.tmap.bam
    else:
        return filename.split('/')[-1].split('_')[0]

bamfiles = sorted(glob.glob(config['tmap_path']), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(bamfiles, parse_sampleID):
    if key.startswith(config['cohort']):   # comment this out if you want to keep blanks
        d[key] = list(value)

# We now have a dictionary with the sample ID as the keys and a list of 
# paths to the bam(s) as the value.
# e.g. d = {sample01: ['a/sample01.bam', 'b/sample01.bam'],
#           sample02: ['a/sample02.bam', 'c/sample02.bam'],
#           sample03: ['d/sample03.bam']}

sampleIDs = d.keys()

# These rules run on the host node and are not submitted to the cluster.
localrules: all

include: 'qc_Snakefile' # creates the fastqc summaries and multiqc report
include: 'annotation_Snakefile' # annotates vcf and creates snpeff multiqc report
include: 'translation_Snakefile' # creates amino acid fasta files
include: 'coverage_Snakefile' # creates a heat map of read coverage


#--------------------------------------------------------------------------
rule all:
    input:
        config['deliver_proj'] + '.fasta',
        'type_summary.tsv',
        'multiqc/fastqc_report.html', # qc_Snakefile
        'multiqc/snpeff_report.html', # annotation_Snakefile
        config['deliver_proj'] + '_aa.fasta', # translation_Snakefile
        'all_ptrim_heatmap.png', # coverage_Snakefile
        'all_raw_heatmap.png' # coverage_Snakefile

def link_rule_input_files(wildcards):
    return d[wildcards.sampleID]

#--------------------------------------------------------------------------
rule link:
    input: link_rule_input_files
    output: 'bams/{sampleID}.bam'
    params:
        bam = 'temp/{sampleID}/{sampleID}.temp.bam',
        sam = 'temp/{sampleID}/{sampleID}.temp.sam'
    run:
        if (len(input) > 1):
            print(os.path.dirname(params.bam))
            shell('mkdir -p %s' %os.path.dirname(params.bam))
            shell('samtools merge {params.bam} {input}')

            # all SM fields in @RG must be identical
            # create a samfile with the fixed header
            sam = pysam.Samfile(params.bam, 'rb')
            header = sam.header
            for i in header['RG']:
                i['SM'] = wildcards.sampleID
            outfile = pysam.Samfile(params.sam, 'wh', header=header)

            # add the reads from the original merged bam
            shell('samtools view {params.bam} >> {params.sam}')

            # convert back to bam
            shell('samtools view -h -b {params.sam} > {output}')
            shell('rm {params.bam} {params.sam}') 
        else:
            shell('cd bams; ln -s {input} {wildcards.sampleID}.bam && touch -h {wildcards.sampleID}.bam')

#--------------------------------------------------------------------------
rule mapq_filter:
    input: rules.link.output
    output: 'mapq_filter/{sampleID}.filtered.bam'
    threads: 2
    params: 
        mapq = int(config["aq_filter"]),
        temp = 'temp/{sampleID}/{sampleID}.aq.temp',
        pre =  'temp/{sampleID}/{sampleID}.sort.temp'
    run:
        shell('mkdir -p %s' %os.path.dirname(params.temp))
        shell('samtools view -h -q {params.mapq} {input} | samtools view -bT {hpv_ref} -o {params.temp}')
        shell('samtools sort -o {output} -@ {threads} -T {params.pre} {params.temp}')
        shell('samtools index {output}; rm {params.temp}')

#--------------------------------------------------------------------------
# Note that TMAP automatically left aligns gaps and indels, so the GATK step
# from the original pipeline was removed.
#--------------------------------------------------------------------------
rule variant_call:
    #input: rules.left_align.output
    input: rules.mapq_filter.output # testing skipping left align
    output: 
        'tvc/{sampleID}/TSVC_variants.vcf',
        'tvc/{sampleID}/{sampleID}.ptrim.bam'
    threads: 2
    params:
        pipe = config["vc_pipe"],
        out = ' tvc/{sampleID}',
        param = config["vc_param"],
        vc_bin = config["vc_bin"],
    run:
        shell('python {params.pipe} \
        --input-bam {input} \
        --postprocessed-bam {output[1]} \
        --primer-trim-bed {amp_bed} \
        --reference-fasta {hpv_ref} \
        --num-threads {threads} \
        --output-dir {params.out} \
        --parameters-file {params.param} \
        --bin-dir {params.vc_bin} \
        --region-bed {len_bed}')


rule adjust_padding:    # this rule doesn't need to run on hpv16 - edit out someday
    input: rules.variant_call.output[0]
    output: 'tvc_vcf/{sampleID}.tvc_no_pad.vcf'
    params: temp = '{sampleID}.temp.vcf'
    run:
        vcf = open(input[0], 'r')
        outfile = open(output[0], 'w')
        need_sort = False

        for line in vcf:
            if line.startswith('#'):
                outfile.write(line)
            else:
                type = line.split()[0]
                hpv_len = bed[type]
                loc = line.split()[1]
                if int(loc) > hpv_len:
                    new_loc = int(loc) - hpv_len
                    outfile.write(line.replace(loc, str(new_loc), 1))
                    need_sort = True
                else:
                    outfile.write(line)
        vcf.close()
        outfile.close()

        if need_sort == True:
            shell('vcf-sort -c {output} > {params.temp}')
            shell('mv {params.temp} {output}')

#--------------------------------------------------------------------------
rule hpv_bam:   # removes human reads
    input: 'tvc/{sampleID}/{sampleID}.ptrim.bam'
    output: 'ptrim_hpv/{sampleID}.hpv.bam'
    run:
        shell('samtools view -h -L {len_bed} {input} | samtools view -bS -o {output}')
        shell('samtools index {output}')

#--------------------------------------------------------------------------
rule pileup:
    input: rules.hpv_bam.output
    output: 'pileup/{sampleID}.pileup'
    run:
        shell('samtools mpileup -f {hpv_ref} -l {len_bed} {input} > {output}')

#--------------------------------------------------------------------------
rule fasta:
    input:
        rules.pileup.output,
        rules.adjust_padding.output
    output:
        'fasta/{sampleID}_HPV%s.fasta' %hpv_type
    run:
        # note vcf header is 0 because it comes after we've skipped 70 lines (as opposed to being the actual line 71)
        df = pandas.read_table(input[1], skiprows=70, header=0)    
        df2 = pandas.read_table(input[0], names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
        types = list(set(df['#CHROM'].tolist() + df2['chrom'].tolist()))

        # create a fasta file for each HPV type found in the sample
        for hpv in types:
            print(hpv)
            num = hpv.replace('HPV', '').replace('_Ref', '')  # convert HPV16_Ref to 16
            # pull out sequence for each HPV type
            # someday make this more efficient by reading them all ahead of time
            fa = config['hpv_reference']['single'] %(num, num)
            fa_handle = open(fa, 'r')

            # only high risk types need fastas
            if os.path.isfile(fa) == False:
                continue

            #record = SeqIO.parse(fa_handle, 'fasta').next()   # just read first record
            # TODO - SeqIOparse().next is throwing an error about no "next" attribute.
            seq = ''
            for record in SeqIO.parse(fa_handle, 'fasta'):
                seq = str(record.seq)
                break # this also takes just the first record, it's just longer
            fa_handle.close()

            # now start looking for SNPs and deletions as per original pipeline
            dt = df[df['#CHROM'] == hpv].copy()
            newseq = seq[:len(seq)-400] # start with the ref sequence and add SNPs
            for idx, row in dt.iterrows():
                (pos, ref, alt) = (int(row['POS']), row['REF'], row['ALT'].split(',')[0])
                # Look for SNPs
                # TODO:  add test for TVC REF matching the REF in the seq string            
                if len(ref) == len(alt):
                    counter = 0
                    while counter < len(ref):
                        newseq = newseq[:pos-1+counter] + alt[counter] + newseq[pos+counter:]
                        counter += 1

                # Look for deletions
                elif (len(ref) > len(alt)) and (len(alt) > 1):
                    counter = 1 # the first nt is the same as the reference, so start at +1
                    while counter < len(ref):
                        newseq = newseq[:pos-1] + '-' + newseq[pos:]
                        counter += 1

                # Skip insertions and other types of variation
                else:
                    continue

            # now check pileup to make sure there was enough coverage at each location
            dp = df2[df2['chrom'] == hpv].copy()

            # Account for zero coverage in pileup (position is missing) by adding zeros
            dp.set_index('pos', inplace=True)
            allpos = list(range(len(seq)+1))
            allpos.pop(0)
            dp = dp.reindex(allpos).fillna(0) 
            dp['cov'] = dp['cov'].astype(int)

            # take padding into account
            dp['adj_pos'] = dp.index
            dp['adj_pos'] = dp['adj_pos'].apply(lambda x: (int(x) - (len(seq)-400)) if x > (len(seq)-400) else int(x))
            x = dp.groupby('adj_pos')['cov'].sum()

            # base calls with less than min_depth are called as N
            x = x[x < config['min_read']]
            for pos, depth in x.iteritems():
                newseq = newseq[:pos-1] + 'N' + newseq[pos:]

            # output a fasta file for each type in each sample
            outfile = open('fasta/%s_HPV%s.fasta' %(wildcards.sampleID, num), 'w')
            outfile.write('>%s_HPV%s\n' %(wildcards.sampleID, num))
            outfile.write(newseq + '\n')
            outfile.close()


#--------------------------------------------------------------------------
rule fasta_cat:
    input: expand('fasta/{sampleID}_HPV%s.fasta' %hpv_type, sampleID=sampleIDs)
    output: config['deliver_proj'] + '.fasta'
    run:
        shell('cat {input} > {output}')

#--------------------------------------------------------------------------
rule type_summary:
    input: expand('pileup/{sampleID}.pileup', sampleID=d.keys())
    output: 'type_summary.tsv'
    run:
        stacks = []
        # iterate through each sample pileup and pull out which types were found
        for sample in input:
            df = pandas.read_table(sample, names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
            df['sampleID'] = sample.split('/')[-1].split('.')[0]
            # note this does not consider padding yet!
            df = df[df['cov'] >= int(config['min_read'])]
            x = df.groupby(['sampleID', 'chrom'])['pos'].count()
            y = x.unstack()
            stacks.append(y) 

        dfs = pandas.concat(stacks).fillna(0)
        dfs.to_csv(output[0], sep='\t')






#--------------------------------------------------------------------------
# convert existing qc_report.sh to snakefile
#onsuccess:
#    # run the qc scripts and copy deliverables to HPV Consortium
#    cov_dev = config["cov_dev"]
#    deliver = config["deliver_proj"]
#    shell('bash qc_report.sh {cov_dev} . {amp_bed} {hg_fasta} {deliver}')
