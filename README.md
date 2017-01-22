# delmh
Find microhomology around deletion

Example:

python delmh/vci2mh.py --vci-type g --vci-s data/all_32_projects/mutect_maf/*.maf --ref-fa hg38.fa --o-func-str "lambda vci_path: os.path.splitext(vci_path)[0].replace('mutect_maf', 'mutect_o') + args.o_suf"