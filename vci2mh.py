#!python

# Xiaoqing Rong-Mullins

from find_mh import * # module search path `sys.path[0]` is the directory where this script is located

import argparse
from Bio import SeqIO
from glob import glob


# print(os.getcwd())


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


# # process

ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_fa, 'fasta'))

for vci_path in vci_list:
    var_table = pd.read_table(vci_path, comment='#')
    dm_out = find_mh_t2t(var_table, ref_dict, vci_cols)
    dm_out.to_csv(o_func(vci_path), index=False)
#


