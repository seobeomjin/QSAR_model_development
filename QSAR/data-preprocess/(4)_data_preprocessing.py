# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 17:49:07 2019

@author: phdgil
"""

import pandas as pd
import numpy as np

wdir ='/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = '5-alpha-reductase-maccs-stdval.csv'

df = pd.read_csv(wdir+csv)
df_colnames = df.columns[2:]

#Remove features with low variance
remove_cols = []
for c in df_colnames:
    if df[c].std() < 0.01:
        remove_cols.append(c)
        
df = df.drop(remove_cols, axis=1)
corr = df.corr().iloc[0]

df.to_csv(wdir+'5-alpha-reductase-maccs-remcols-stdval.csv', index=False)