#!/usr/bin/env python3

import os
import pandas
import re

#############
# FUNCTIONS #
#############

def fix_sample_names(sample_name):
    if sample_name.endswith('*'):
        return(re.sub('^R(\d+).*', 'Rpoa\\1', sample_name))
    elif sample_name.endswith('-H'):
        return(re.sub('-H$', '', sample_name))
    else:
        return(sample_name)


def find_key_files(wildcards):
    # get a list of key files
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    all_key_files = []
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.txt'):
                all_key_files.append(os.path.join(dirpath, filename))
    return(all_key_files)


def find_read_file(fc, data_dir):
    # get a list of key files
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz') and fc in filename:
                return(os.path.join(dirpath, filename))

def list_fc_names(data_dir):
    # get a list of key files
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    all_fcs = []
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.txt'):
                all_fcs.append(re.sub('.txt', '', filename))
    return(all_fcs)


def read_key_and_write_config(key_file, outdir):
    # generate filename
    my_filename = re.sub('\\.txt$', '_config', os.path.basename(key_file))
    my_path = os.path.join(outdir, my_filename)
    # read keyfile
    key_data = pandas.read_csv(key_file, delimiter='\t')
    key_data.dropna(how='all', inplace=True)
    # remove Ophir samples
    filtered_keydata = key_data[
        list(not x.startswith('O') for x in key_data['sample'])]
    # remove asterisks and '-H' from sample names
    filtered_keydata['sample_name'] = filtered_keydata['sample'].apply(
        fix_sample_names)
    # write the two columns of interest to an output file
    subset = filtered_keydata[['barcode', 'sample_name']]
    subset.to_csv(my_path,
                  sep='\t',
                  header=False,
                  index=False)




###########
# GLOBALS #
###########

data_dir = 'data/asw_para_matched'

#########
# RULES #
#########

rule target:
    input:
        dynamic(expand('output/020_demux/{fc}/{{individual}}.fq.gz',
                fc = list_fc_names(data_dir)))

for fc in list_fc_names(data_dir):
    my_config_file = 'output/010_config/{}_config'.format(fc)
    my_read_file = find_read_file(fc, data_dir)
    my_outdir = 'output/020_demux/{}'.format(fc)
    my_output_files = 'output/020_demux/{}/{{individual}}.fq.gz'.format(fc)
    rule:
        input:
            config = my_config_file,
            reads = my_read_file
        output:
            dynamic(my_output_files)
        params:
            outdir = my_outdir
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


rule generate_config_files:
    input:
        find_key_files
    threads:
            1
    output:
        'output/010_config/{fc_name}_config'
    params:
        outdir = 'output/010_config'
    run:
        for key_file in input:
            read_key_and_write_config(key_file, params.outdir)
            