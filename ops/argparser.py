#
# Copyright (C) 2020 Xiao Wang
# Email:xiaowang20140001@gmail.com wang3702@purdue.edu
#
import argparse
import json
from collections import OrderedDict
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',type=int,required=True,help='control mode')
    parser.add_argument('-F',type=str, help='input map path')#File path for decoy dir
    parser.add_argument("-P",type=str,help="directory or zipped file of Single-Chain PDB files")
    parser.add_argument("-M",type=str,help="txt file path which records protein information")
    parser.add_argument("--resolution",type=str,help="specify the resolution to use different model")
    parser.add_argument("--config",type=str,default=None,help="specifying the config path")
    parser.add_argument("--gpu",type=str,default=None,help="specify the gpu we will use")
    parser.add_argument("--output",type=str,help="Output directory")
    parser.add_argument("--contour",type=float,default=0,help="Contour level for input map, suggested 0.5*[author_contour]. (Float), Default value: 0.0")
    args = parser.parse_args()
    # remove comments starting with '//'
    json_str = ''
    opt_path = args.config
    params = vars(args)
    if opt_path is not None:
        with open(opt_path, 'r') as f:
            for line in f:
                line = line.split('//')[0] + '\n'
                json_str += line
        opt = json.loads(json_str, object_pairs_hook=OrderedDict)

        for key in opt:
            params[key]=opt[key]
    return params
