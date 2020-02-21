# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 16:28:26 2019

@author: phdgil
"""

from rdkit.Chem import MACCSkeys                       #MACCS Keys
from rdkit import Chem, DataStructs
import numpy as np
import pandas as pd

def fingerprint_bitgenerator(fps):
    np_fps = []
    
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp,arr)
        np_fps.append(arr)
        
    return np_fps

wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = 'final_ic50.csv'

df = pd.read_csv(wdir+csv)
smi = df['Canonical Smiles']
sd = [Chem.MolFromSmiles(m) for m in smi]

maccfps = [MACCSkeys.GenMACCSKeys(m) for m in sd]
maccfps_bit = fingerprint_bitgenerator(maccfps)
fp_df = pd.DataFrame(maccfps_bit,columns=list(range(len(maccfps_bit[0]))))

df_id = df[['Molecule','Standard Value']]

#Data table preparation
df_label_fp = pd.concat([df_id, fp_df], axis=1)
df_label_fp.to_csv(wdir+'5-alpha-reductase-maccs-stdval.csv',index=False)

#PCA: Data plotting
'''
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca = PCA(n_components = 2).fit(fp_df)
fp_pca = pca.transform(fp_df)
fp_pca_df = pd.DataFrame(fp_pca, columns=['pca1','pca2'])
df_label_pca = pd.concat([df_id, fp_pca_df], axis=1)

inhib = df_label_pca[df_label_pca['label']==1]
non = df_label_pca[df_label_pca['label']==0]

plt.scatter(inhib['pca1'],inhib['pca2'],c='red',label='inhibitor')
plt.scatter(non['pca1'],non['pca2'],c='blue',label='non-inhibitor')
plt.title('5-alpha-reductase inhibition')
plt.xlabel('MACCS pca axis 1')
plt.ylabel('MACCS pca axis 2')
plt.legend(loc=4)
plt.savefig(wdir+'visual/'+'Data_distribution_MACCS_pca_stdval.png')
'''