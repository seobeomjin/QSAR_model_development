# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:37:12 2019

@author: phdgil
"""
# RandomForestRegressor

import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score

#data load
wdir = '/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
csv = '5-alpha-reductase-maccs-remcols-stdval.csv'

df = pd.read_csv(wdir+csv)
df = df.dropna(axis=0)
y = df['Standard Value']
x = df.drop(['Molecule', 'Standard Value'], axis=1)

#data split
x_train, x_test, y_train, y_test = train_test_split(x,y,test_size = 0.1)
train_size = len(x_train)

#model fit
RFregr = RandomForestRegressor()
RFregr.fit(x_train,y_train)
cv = 5

#metric
mse = cross_val_score(RFregr, 
                        x_train, 
                        y_train, 
                        scoring='neg_mean_squared_error',  
                        cv=cv )
mse *= -1
rmse  = np.sqrt(mse)

print('feature_importances',RFregr.feature_importances_)
print('RandomForestRegression RMSE : ', rmse)
print('RandomForestRegression RMSE.mean : ' , np.sum(rmse)/cv)

