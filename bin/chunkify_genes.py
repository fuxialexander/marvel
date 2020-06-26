#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--promoter-gene-pair-file", default="promoter_gene_pair.txt",
                    required=True, type=str, help="output directory")
parser.add_argument('-n', "--chunk-size", default=100, type=int, required=True,
                    help="Chunk size of gene splits")
parser.add_argument('-o', "--output-path", default=".",
                    required=True, type=str, help="output directory")
args = parser.parse_args()


pair = pd.read_csv(args.promoter_gene_pair_file, '\t',
                   header=None, names=['promoter', 'gene'])


chunk_list = [pd.DataFrame(x) for i, x in pair.groupby('gene', as_index=False)]


chunk_index = np.arange(len(chunk_list))


chunk_index_list = np.array_split(
    chunk_index, max(2, len(chunk_index)//args.chunk_size))


for i, g in enumerate(chunk_index_list):
    fname = args.output_path+'/promoter_gene_pair_' + str(i) + '.txt'
    gi = pd.concat([chunk_list[j] for j in g], 0)
    gi.to_csv(fname, sep='\t', header=False, index=False)
