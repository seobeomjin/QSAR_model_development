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

    mse = cross_val_score(rfe, 
                        x_train, 
                        y_train, 
                        scoring='neg_mean_squared_error',  
                        cv=cv )
    mse *= -1
    rmse  = np.sqrt(mse)

    print('{} SVR with {} features RMSE : '.format(kernel, f), rmse)
    print('{} SVR with {} features RMSE.mean : '.format(kernel, f) , np.sum(rmse)/cv)
    


