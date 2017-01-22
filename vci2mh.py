#!python

# Xiaoqing Rong-Mullins

import argparse
import numpy as np
import os
import pandas as pd

from Bio import SeqIO
from glob import glob


# # Parser

parser = argparse.ArgumentParser(description='Find microhomology.')

# variant call input (vci)
parser.add_argument(
    '--vci-type', default='g',
    help='\
        Type of information specified by --vci-s. Available options: \
            \'d\' (directory in which all the items are vci files), \
            \'f\' (filename of vci), \
            \'g\' or anything else (`pathname` parameter for `glob` to find vci files). \
            Default: \'g\'.\
    '
)
parser.add_argument(
    '--vci-s', action='append',
    help='\
        Information about variant call input (vci) paths. Must be specified, no default. \
        Multiple vci paths can be specified by "--vci-s path1 --vci-s path2 ...".\
    '
)
parser.add_argument(
    '--vci-cols-str',
    default='\
Hugo_Symbol,Entrez_Gene_Id,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,Variant_Type,\
Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Tumor_Sample_Barcode,Matched_Norm_Sample_Barcode\
',
    help='Headers of columns from the vci files to be included in output, separated by --vci-cols-sep. Default has 13 headers.'
)
parser.add_argument(
    '--vci-cols-sep', default=',',
    help='Delimiter to separate column headers in --vci-cols-str using `str.split` method. Default: , (comma)'
)

# reference genome .fasta
parser.add_argument('--ref-fa', help='Path to reference genome .fasta. Must be specified, no default.')

# output
parser.add_argument(
    '--o-path',
    help='Path to one output file, valid when there is only one vci file; if not specified, use `o_func` get it'
)
parser.add_argument('--o-pref', default='', help='Prefix for output path. Default is an empty string.')
parser.add_argument('--o-suf', default='.csv', help='Suffix for output path. Default: \'.csv\'.')
parser.add_argument(
    '--o-func-str', 
    default='lambda vci_path: args.o_pref + os.path.splitext(vci_path)[0] + args.o_suf',
    help='\
        A string to be converted by `eval` into a function to obtain `output_path`. \
        Default is to concatenate --o-pref, the root of vci path (split by `os.path.split`) and --o-suf together.\
    '
)


# # Parse arguments

args = parser.parse_args()

# vci
if args.vci_type == 'f':
    vci_list = args.vci_s
else: # if args.vci_type == 'd', every item in each `dir_` will be used as an vci file
    vci_list = []
    for item in args.vci_s:
        glob_pathname = item + '/*' if args.vci_type == 'd' else item
        vci_list += glob(glob_pathname)
vci_cols = args.vci_cols_str.split(args.vci_cols_sep)

# o_func
print(args.o_func_str)
if args.o_path != None and len(vci_list) == 1:
    o_func = lambda x: args.o_path
else:
    o_func = eval(args.o_func_str)
#


# # functions

coord2pyidx = -1

def find_mh_1d(direction, seq, start_coord, end_coord): # python_index = coord - 1; coord = python_index + 1
    
    if direction == 'up':
        step = -1
        offset = 0
    elif direction == 'down':
        step = 1
        offset = 1
    else:
        raise Exception('`direction` must be either \'up\' or \'down\'.')
    
    mh = ''
    left_idx = init_left_idx = int(start_coord) + coord2pyidx - 1 + offset # numpy.int64 is "Invalid index"
    right_idx = init_right_idx = int(end_coord) + coord2pyidx + offset
    while seq[left_idx] == seq[right_idx]:
        mh += seq[left_idx]
        left_idx += step
        right_idx += step
    if left_idx == init_left_idx:
        found = False
        left_coords = right_coords = []
    else:
        found = True
        coords = [idx - coord2pyidx for idx in [init_left_idx, left_idx - step, init_right_idx, right_idx - step]]
        left_coords = coords[:2]
        right_coords = coords[2:]
        if direction == 'up':
            mh = mh[::-1]
            left_coords = left_coords[::-1]
            right_coords = right_coords[::-1]
    
    return {'found': found, 'mh': mh, 'left_coords': left_coords, 'right_coords': right_coords}
        
    
def find_mh_2d(seq, start_coord, end_coord, out_type='dict'): # idx is python style
    
    directions = ['up', 'down']
    
    mh_dict = {direc: find_mh_1d(direc, seq, start_coord, end_coord) for direc in directions}
    
    mh_list = []
    for direc in directions:
        mh_dict_direc = mh_dict[direc]
        mh_srs_dict_direc = {key: mh_dict_direc[key] for key in ['found', 'mh']}
        for position in ['left', 'right']:
            coords = mh_dict_direc[position + '_coords']
            if len(coords) == 0:
                coords = [np.nan] * 2
            pos_dict = {position + '_start': coords[0], position + '_end': coords[1]}
            mh_srs_dict_direc.update(pos_dict)
        mh_srs_direc = pd.Series(mh_srs_dict_direc)
        mh_srs_direc.index = ['_'.join([direc, idx]) for idx in mh_srs_direc.index]
        mh_list.append(mh_srs_direc)
    mh_srs = pd.concat(mh_list, axis=0)
    
    if out_type == 'dict':
        return mh_dict
    elif out_type == 'srs':
        return mh_srs
    else:
        raise Exception('`out_type` must be either \'dict\' or \'srs\'.')
#


# # process

ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_fa, 'fasta'))

for vci_path in vci_list:
    var_table = pd.read_table(vci_path, comment='#')
    dels = var_table.loc[var_table['Variant_Type'].apply(lambda val: val == 'DEL')]
    del_mh = dels.apply(lambda row: find_mh_2d(
        ref_dict[row['Chromosome']], row['Start_Position'], row['End_Position'], out_type='srs'
    ), axis=1)
    del_out = pd.concat([dels[vci_cols], del_mh], axis=1)
    del_out.to_csv(o_func(vci_path), index=False)
#


