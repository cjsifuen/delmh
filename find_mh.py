# find_mh

import numpy as np
import os
import pandas as pd


def find_mh_1d(direction, seq, start_coord, end_coord, coord2pyidx = -1): # python_index = coord - 1; coord = python_index + 1
    
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
    while (
            left_idx >= 0 and right_idx < len(seq) # search within the bounds of reference sequence
        and 
            abs(left_idx - init_left_idx) <= end_coord - start_coord # one member of homology is within the deletion region
        and 
            seq[left_idx] == seq[right_idx] # homology detected
        ):
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

vci_cols_default = ['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type']

def find_mh_t2t(var_table, ref_dict, vci_cols=vci_cols_default): # from input table to output table
    dels = var_table.loc[var_table['Variant_Type'].apply(lambda val: 'DEL' in val)]
    del_mh = dels.apply(lambda row: find_mh_2d(
        ref_dict[row['Chromosome']], row['Start_Position'], row['End_Position'], out_type='srs'
    ), axis=1)
    return pd.concat([dels[vci_cols], del_mh], axis=1)
#




