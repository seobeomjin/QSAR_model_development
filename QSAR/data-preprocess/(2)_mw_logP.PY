# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 16:11:32 2019

@author: phdgil
"""

import matplotlib.pyplot as plt
import pandas as pd

wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = 'final_ic50.csv'

df = pd.read_csv(wdir+csv)

inhib = df[df['label']==1]
non = df[df['label']==0]

plt.scatter(inhib['ALogP'],inhib['Molecular Weight'],c='red',label='inhibitor')
plt.scatter(non['ALogP'],non['Molecular Weight'],c='blue',label='non-inhibitor')
plt.title('5-alpha-reductase inhibition')
plt.xlabel('logP')
plt.ylabel('MW')
plt.legend(loc=4)
plt.savefig(wdir+'Data_distribution_logP_MW.png')