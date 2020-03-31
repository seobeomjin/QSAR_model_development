# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:37:12 2019

@author: phdgil
"""
# SVR , cross-validation , MSE 

import random
import pandas as pd
from sklearn.svm import SVR
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score
import numpy as np 

#data load
wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = '5-alpha-reductase-maccs-remcols-stdval.csv'

df = pd.read_csv(wdir+csv)
df = df.dropna(axis=0)
y = df['Standard Value']
x = df.drop(['Molecule', 'Standard Value'], axis=1)

#hyperparam
test_size = 0.1
kernel = 'poly' #if you want, you can try other kernels like 'linear','rbf'
cv = 5

#data split
x_train, x_test, y_train, y_test = train_test_split(x,y,test_size = test_size)

#model fit
svr = SVR(kernel=kernel, gamma='auto')
svr.fit(x_train,y_train)

#metric
mse = cross_val_score(svr, 
                        x_train, 
                        y_train, 
                        scoring='neg_mean_squared_error',  
                        cv=cv )
mse *= -1
rmse  = np.sqrt(mse)

print('{} SVR RMSE : '.format(kernel), rmse)
print('{} SVR RMSE.mean : '.format(kernel) , np.sum(rmse)/cv)
