# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:37:12 2019

@author: phdgil
"""
import random
import pandas as pd
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt

import numpy as np 

def SquredError (x,y):
    se = ((x-y)*(x-y))
    return se

def bootstrapping (datasize,x, y):
    #Random sampling with duplication
    train = []
    for i in range(datasize):
        train.append(random.randrange(datasize))
    
    #Find indices not sampled
    test = []
    for i in range(df.shape[0]):
        if i not in train:
            test.append(i)
    
    #Subsetting dataframe from sampled and unsampled data
    train_x = x.iloc[train]
    train_y = y.iloc[train]
    
    test_x = x.iloc[test]
    test_y = y.iloc[test]
    
    return train_x, train_y, test_x, test_y

wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = '5-alpha-reductase-maccs-remcols.csv'

df = pd.read_csv(wdir+csv)
y = df['label']
x = df.drop(['Molecule', 'label'], axis=1)
datasize = df.shape[0]
sampling = 50

svc = SVC(kernel='linear', gamma='auto')
feature = range(10,30,2)  
acc = []
val = []
for f in feature:
    print ('feature :', f)
    rfe = RFE(svc, n_features_to_select=f, step=1)  
    
    #Validation with bootstrapping
    boot_acc = []
    boot_val = []
    for i in range(sampling):
        train_x, train_y, test_x, test_y = bootstrapping(datasize, x, y)
        rfe.fit(train_x,train_y) 
        boot_acc.append(rfe.score(train_x, train_y)) 
        boot_val.append(rfe.score(test_x, test_y))
         
    acc.append(sum(boot_acc)/len(boot_acc))  
    val.append(sum(boot_val)/len(boot_val))


plt.plot(feature, acc, 'r-', label='train')
plt.plot(feature,val, 'b-', label='test')
plt.legend(loc=4)
plt.xlabel('Number of features')
plt.ylabel('Accuracy')
plt.title('SVC model with linear kernel (10~30 features)')
plt.savefig(wdir+'svc_feature_selection_10-30.png')


