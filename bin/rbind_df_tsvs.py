#!/usr/bin/env python

import sys
import pandas as pd

out_tsv = sys.argv[1]
in_tsvs = sys.argv[2:]

df = pd.concat([pd.read_csv(x, sep='\t') for x in in_tsvs],
               axis=0, ignore_index=True)

df.to_csv(
    out_tsv,
    sep='\t', index=False
)
