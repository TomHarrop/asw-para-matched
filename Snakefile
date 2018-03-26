#!/usr/bin/env python3

import os
import pandas
import re

#############
# FUNCTIONS #
#############

def fix_rpoa_sample_names(sample_name):
    if sample_name.endswith('*'):
        return(re.sub('^R(\d+).*', 'Rpoa\\1', sample_name))
    elif sample_name.endswith('-H'):
        return(re.sub('-H$', '', sample_name))
    else:
        return(sample_name)


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
        fix_rpoa_sample_names)
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
# SETUP #
#########



#########
# RULES #
#########

