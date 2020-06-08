"""
Written by Alexandre Pellan Cheng 2020
WGBS Alignmnent pipeline.

Works for preps done via Swift, Meyer, SRSLY.

This file is the motherboard of the operation. All softwares,
paths, variables and functions are defined here.

rule all is declared here as well.

Used "include" to call all necessary functions to go from zero to hero.
By that, I mean that references are downloaded and indexed and so you only need this file +
a few auxiliary files to run the entire pipeline.

"""
from natsort import natsorted
import itertools
from snakemake.utils import validate, min_version

min_version("5.1.2")
################################################################################
# CONFIG FILE
################################################################################
configfile: 'config.yaml'
################################################################################
# FUNCTIONS
################################################################################
def string_split(string, delimiter, number_of_splits):
    s=string.split(delimiter, number_of_splits)
    return(s)

tissues = natsorted(set(x.split('_')[0] for x in config['REFERENCE_METHYLOMES']))
comparing_groups=list(itertools.combinations(tissues, 2))
comparing_groups=[x[0]+'_'+x[1] for x in comparing_groups]


def get_fastq_reads(wcrds):
    """
    We get the FQs based on if we got SE or PE
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            SAM, PREP, SEQ, CTL_type = line.strip().split('\t')
            if str(wcrds.sample) == SAM:
                seq = SEQ.split('x')[0]
                prep = PREP
    if seq == '1': #single-end
        return[CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz', CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz'] #second R1 is a dumb placeholder
    if seq == '2': #paired-end
        return[CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz', CLUSTER+config['DATA']+'samples/{sample}_R2.fastq.gz']

    raise(ValueError('Could not figure paired end or single end'))

def get_seq_type(wcrds):
    """
    We get the FQs based on if we got SE or PE
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            SAM, PREP, SEQ, CTL_type = line.strip().split('\t')
            if str(wcrds.sample) == SAM:
                seq =SEQ
                prep = PREP
    return[prep, seq]

def get_adapter_file(wcrds):
    """
    For trimming
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            SAM, PREP, SEQ, CTL_type = line.strip().split('\t')
            if str(wcrds.sample) == SAM:
                prep = PREP
    if prep == "MEYER_SSLP":
        adapt = MEYER_ADAPTORS
    if prep == "SRSLY_SSLP":
        adapt = SRSLY_ADAPTORS
    if prep == "SWIFT_ACCEL":
        adapt = SA_ADAPTORS
    return(adapt)



################################################################################
# CLUSTER USED
################################################################################
CLUSTER = config['CLUSTER']
################################################################################
# File locations
################################################################################
SEQUENCING_PREP_TABLE = CLUSTER + config['SEQUENCING_PREP_TABLE']
MEYER_ADAPTORS = CLUSTER + config['MEYER_ADAPTORS']
SRSLY_ADAPTORS = CLUSTER + config['SRSLY_ADAPTORS']
SA_ADAPTORS = CLUSTER + config['SA_ADAPTORS']
HG19METH = CLUSTER + config['HG19METH']
CTLMETH = CLUSTER + config['CTLMETH']
METHYLOME_TABLE = CLUSTER + config['METHYLOME_TABLE']
################################################################################
# Software paths
################################################################################
SPRING = CLUSTER + config['SPRING']
BISMARKINDEX = CLUSTER + config['BISMARKINDEX']
BBDUK = CLUSTER + config['BBDUK']
SGREP = CLUSTER + config['SGREP']
SEQTK = CLUSTER + config['SEQTK']
HMMCOPY_FA_BW = CLUSTER + config['HMMCOPY_FA_BW']
HMMCOPY_MAP = CLUSTER + config['HMMCOPY_MAP']
HMMCOPY_GC = CLUSTER + config['HMMCOPY_GC']
HMMCOPY_READCOUNTER = CLUSTER + config['HMMCOPY_READCOUNTER']
BISMARK = CLUSTER + config['BISMARK']
METHEXT = CLUSTER + config['METHEXT']
RMDUPS = CLUSTER + config['RMDUPS']
REPAIR = CLUSTER + config['REPAIR']
METHPIPETOMR = CLUSTER + config['METHPIPETOMR']
METHPIPEBSRATE = CLUSTER + config['METHPIPEBSRATE']
HSBLASTN = CLUSTER + config['HSBLASTN']
GRAMMY_GDT = CLUSTER + config['GRAMMY_GDT']
GRAMMY_RDT = CLUSTER + config['GRAMMY_RDT']
GRAMMY_PRE = CLUSTER + config['GRAMMY_PRE']
GRAMMY_EM = CLUSTER + config['GRAMMY_EM']
GRAMMY_POST = CLUSTER + config['GRAMMY_POST']
GRAMMY_REF_FASTA = CLUSTER + config['GRAMMY_REF_FASTA']
bigWigToBedGraph = config['bigWigToBedGraph']
CROSSMAP = CLUSTER + config['CROSSMAP']
METILENE = CLUSTER + config['METILENE']
TRIM_GALORE = CLUSTER + config['TRIM_GALORE']
################################################################################
# Variables
################################################################################
AUTOSOMALCHROMO = config['AUTOSOMALCHROMO']
KNOWNCHROMO = config['KNOWNCHROMO']
################################################################################
# Threads
################################################################################
indexing_threads=15
fastqc_threads=1
trim_threads=4
alignment_threads=4
filter_bam_threads=4
methylation_extraction_threads=4
blast_threads=8
################################################################################
# RULES
################################################################################
rule all:
    input:
        expand('sample_output/fastqc/{sample}_{read}_fastqc.html', sample = config['SAMPLES'], read = ['R1', 'R2']),
        expand('sample_output/stats/{sample}_mapping_stats.txt', sample=config['SAMPLES']),
        expand('sample_output/Lengths/{sams_ctls}_FragsHistogram.pdf', sams_ctls = config['SAMPLES']),
        expand('sample_output/conversion_rates/{sample}.bsconversion', sample = config['SAMPLES']),
        #expand('sample_output/donor_fraction/{samples}.XY', samples = config['SAMPLES']),
        expand('sample_output/mbias/{sample}.pdf', sample = config['SAMPLES']),
        expand('sample_output/classic_tissues_of_origin{dir}/{sample}.tsv', dir = ['_trim', '_untrim'], sample = config['SAMPLES']),
        expand('sample_output/total_abundance/{sample}.txt', sample = config['SAMPLES']),
        #expand('sample_output/grammy/{sams_ctls}/{sams_ctls}.grammy.tab', sams_ctls = config['SAMPLES'], conversion=['CT', 'GA'])

#Rules that create references (generally only need to be run once)
#include: 'rules/references/get_and_index_hg19.smk'
#include: 'rules/references/index_mic_ctl.smk'
#include: 'rules/references/get_and_index_phix.smk'
#include: 'rules/references/get_and_index_GRAMMy.smk'
#include: 'rules/references/prepare_donorfraction.smk'

#include: 'rules/reference_methylomes/download_methylomes.smk'
#include: 'rules/reference_methylomes/format_methylomes.smk'
#include: 'rules/reference_methylomes/create_methylmatrix.smk'
#include: 'rules/reference_methylomes/get_DMRs2.smk'



include: 'rules/qc/fastqc.smk'
include: 'rules/alignment/trim.smk'
include: 'rules/alignment/alignment.smk'
include: 'rules/stats/mapping_stats.smk'
include: 'rules/stats/total_abundance.smk'
include: 'rules/alignment/length_profiles.smk'
include: 'rules/alignment/donor_fraction.smk'
include: 'rules/methylation/estimate_BSconversion.smk'
include: 'rules/methylation/methylation_extraction.smk'
include: 'rules/methylation/tissues_of_origin.smk'
include: 'rules/metagenome/phix_alignment.smk'
include: 'rules/metagenome/decontaminate.smk'
include: 'rules/metagenome/hsblastn.smk'
include: 'rules/metagenome/grammy.smk'
