#!/usr/bin/env python3

import os
import pandas

#############
# FUNCTIONS #
#############


def all_read_files(data_dir):
    '''
    Return the fastq.gz files that matches wildcards.fc
    '''
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    all_read_files = []
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                all_read_files.append(os.path.join(dirpath, filename))
    return(all_read_files)


def generate_fc_dicts(key_file, data_dir):
    '''
    Parse the key_file and return dicts
    '''
    # initialise the dictionaries
    fc_to_indiv = {}
    fc_to_readfile = {}
    # find files
    read_files = all_read_files(data_dir)
    # read the key data
    key_data = pandas.read_csv(key_file)
    grouped_key_data = key_data.groupby('key')
    # populate dicts
    for name, group in grouped_key_data:
        fc_to_indiv[name] = sorted(set(x for x in group['sample_name']))
        fc_to_readfile[name] = [x for x in read_files
                                if name in os.path.basename(x)][0]
    return(fc_to_indiv, fc_to_readfile)


def read_keydata_and_write_config(key_file, outdir):
    '''
    Parse the keyfile and generate config files in outdir
    '''
    # read the key data
    key_data = pandas.read_csv(key_file)
    grouped_key_data = key_data.groupby('key')
    for name, group in grouped_key_data:
            config_file = os.path.join(outdir, '{}_config'.format(name))
            subset = group[['barcode', 'sample_name']]
            if len(subset) > 0:
                subset.to_csv(config_file,
                              sep='\t',
                              header=False,
                              index=False)


###########
# GLOBALS #
###########

data_dir = 'data/asw_para_matched'
key_file = 'data/asw_para_matched/combined_key_data.csv'

#########
# SETUP #
#########

fc_to_indiv, fc_to_readfile = generate_fc_dicts(
    key_file,
    data_dir)
all_fcs = list(set(fc_to_indiv.keys()))
all_indivs = sorted(set(y for x in all_fcs for y in fc_to_indiv[x]))

#########
# RULES #
#########

rule target:
    input:
        'output/030_optim/compare_defaults/optimised_samplestats_combined.csv'

rule compare_defaults:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/compare_defaults/optimised_samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/020_demux'
    log:
        'output/logs/030_optim/compare_defaults.log'
    shell:
        'stacks_parameters '
        '--mode compare_defaults '
        '-m 3 '
        '-M 3 '
        '-n 3 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 2 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_n:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/stats_n/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/020_demux'
    log:
        'output/logs/030_optim/optim_n.log'
    shell:
        'stacks_parameters '
        '--mode optim_n '
        '-m 3 '
        '-M 3 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 2 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

rule optim_mM:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/stats_Mm/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/020_demux'
    log:
        'output/logs/030_optim/optim_mM.log'
    shell:
        'stacks_parameters '
        '--mode optim_Mm '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 2 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_setup:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/filtering/replicate_1_popmap.txt'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/020_demux'
    log:
        'output/logs/030_optim/optim_setup.log'
    shell:
        'stacks_parameters '
        '--mode setup '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

rule generate_popmap:
    input:
        key_file = key_file
    output:
        popmap = 'output/010_config/full_popmap.txt'
    threads:
        1
    run:
        key_data = pandas.read_csv(input.key_file)
        subset = key_data.loc[key_data['population'] != 'GBSNEG',
                              ['sample_name', 'population']]
        subset.to_csv(output.popmap,
                      sep='\t',
                      header=False,
                      index=False)


for fc in all_fcs:
    rule:
        input:
            config = 'output/010_config/{}_config'.format(fc),
            reads = fc_to_readfile[fc]
        output:
            expand('output/020_demux/{individual}.fq.gz',
                   individual=fc_to_indiv[fc])
        params:
            outdir = 'output/020_demux'
        threads:
            1
        log:
            'output/logs/020_demux/{}.log'.format(fc)
        shell:
            'process_radtags '
            '-f {input.reads} '
            '-i gzfastq -y gzfastq '
            '-b {input.config} '
            '-o {params.outdir} '
            '-c -q '
            '-t 91 '
            '--inline_null '
            '--renz_1 apeKI --renz_2 mspI '
            '&> {log} '

# 010 generate stacks config
rule generate_config_files:
    input:
        key_file = key_file
    threads:
            1
    output:
        expand('output/010_config/{fc_name}_config',
               fc_name=all_fcs)
    params:
        outdir = 'output/010_config'
    run:
        read_keydata_and_write_config(input.key_file, params.outdir)

