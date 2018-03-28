#!/usr/bin/env python3

import os
import pandas
import re


#############
# FUNCTIONS #
#############

def fix_sample_names(sample_name):
    '''
    Remove -H suffix and/or asterisk from sample_name
    '''
    if sample_name.endswith('*'):
        return(re.sub('^R(\d+).*', 'Rpoa\\1', sample_name))
    elif sample_name.endswith('-H'):
        return(re.sub('-H$', '', sample_name))
    else:
        return(sample_name)


def find_key_files(wildcards):
    '''
    Return a list of all .txt from data_dir
    '''
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    all_key_files = []
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.txt'):
                all_key_files.append(os.path.join(dirpath, filename))
    return(all_key_files)


def find_read_file(fc):
    '''
    Return the fastq.gz file that matches wildcards.fc
    '''
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if (filename.endswith('.fastq.gz')
                    and fc in filename):
                return(os.path.join(dirpath, filename))


def list_fc_names(data_dir):
    '''
    Return a list of  in data_dir without .txt
    '''
    return(
        [os.path.splitext(os.path.basename(x))[0] for x in find_key_files('')])


def read_key_and_write_keydata(key_file, outdir):
    '''
    Read the key_file, fix the sample names, and write a csv of keydata to
    outdir
    '''
    # generate filename
    my_filename = re.sub('\\.txt$', '_keydata.csv', os.path.basename(key_file))
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
    # write the columns of interest to an output file
    keep_columns = ['flowcell', 'lane', 'barcode', 'sample', 'sample_name']
    subset = filtered_keydata[keep_columns]
    output_df = subset.rename(columns={'sample': 'agr_sample_name'})
    output_df.to_csv(my_path,
                     sep=',',
                     header=True,
                     index=False)


def read_keydata_and_write_config(keydata, outdir):
    my_filename = re.sub('_keydata.csv', '_config', os.path.basename(keydata))
    my_path = os.path.join(outdir, my_filename)
    # read keyfile
    key_data = pandas.read_csv(keydata, delimiter=',')
    key_data.dropna(how='all', inplace=True)
    # write the columns of interest to an output file
    keep_columns = ['barcode', 'sample_name']
    output_df = key_data[keep_columns]
    output_df.to_csv(my_path,
                     sep='\t',
                     header=False,
                     index=False)


###########
# GLOBALS #
###########

data_dir = 'data/test_data'

#########
# RULES #
#########

rule target:
    input:
        dynamic(expand('output/020_demux/{fc_name}/{{individual}}.fq.gz',
                fc_name=list_fc_names(data_dir)))
    

# 020 demux
for fc in list_fc_names(data_dir):
    rule:
        input:
            config = 'output/010_config/{}_config'.format(fc),
            reads = find_read_file(fc)
        output:
            dynamic('output/020_demux/{}/{{individual}}.fq.gz'.format(fc))
        params:
            outdir = 'output/020_demux/{}'.format(fc)
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
        key_file = 'output/010_config/{fc_name}_keydata.csv'
    threads:
            1
    output:
        'output/010_config/{fc_name}_config'
    params:
        outdir = 'output/010_config'
    run:
        read_keydata_and_write_config(input.key_file, params.outdir)

rule read_key_data:
    input:
        key_file = 'data/asw_para_matched/{fc_name}.txt'
    threads:
            1
    output:
        'output/010_config/{fc_name}_keydata.csv'
    params:
        outdir = 'output/010_config'
    run:
        read_key_and_write_keydata(input.key_file, params.outdir)

