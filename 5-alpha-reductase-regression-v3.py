# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:37:12 2019

@author: phdgil
"""
# SVR , cross-validation , MSE 
# with Recursive Feature Elimination (RFE) to get the informative features 

import random
import pandas as pd
from sklearn.svm import SVR
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score

import numpy as np 
'''
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
'''

wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = '5-alpha-reductase-maccs-remcols-stdval.csv'

df = pd.read_csv(wdir+csv)
df = df.dropna(axis=0)
y = df['Standard Value']
x = df.drop(['Molecule', 'Standard Value'], axis=1)
#datasize = df.shape[0]
#sampling = 200

#hyperparam for cross_val
test_size = 0.1  
cv = 5  
x_train, x_test, y_train, y_test = train_test_split(x,y,test_size = test_size)  #data split 

#hyperparam for SVR
kernel = 'linear'
svr = SVR(kernel=kernel, gamma='auto')
feature = range(10,30,2)  ### N of features is selected as 10, 12, 14 ,,, to 30 

for f in feature:
    
    print ('feature :', f)
    rfe = RFE(svr, n_features_to_select=f, step=1)  
    ### when N of feature are 10, 12, 14 ,, 30 
    ### take feature selection instantiate to 'rfe' with svc as param
    
    rfe = rfe.fit(x_train,y_train)
    print('N of features : {}\n feature support : {}\n feature ranking : {} '.format(f,rfe.support_,rfe.ranking_))
    


