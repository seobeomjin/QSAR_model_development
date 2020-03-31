from rdkit.Chem import MACCSkeys                       #MACCS Keys
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

class std_tendency():
    
    # if i instantiate something here, 
    # a class variable is shared by all class instances
    
    # >> it would be a 'class attribute'
    
    def __init__(self,wdir,csv,chem_num): 
        # >> it would be a instance attribute 
        self.wdir = wdir
        self.csv = csv
        self.chem_num = chem_num
        
        self.c0, self.c1, self.c2, self.c3, self.c4 = [], [], [], [], []
        self.c0_x, self.c1_x, self.c2_x, self.c3_x, self.c4_x = [], [], [], [], []
        self.c4_stdval, self.c4_stdval_x = [], []
        
        self.df = pd.read_csv(self.wdir+self.csv)
        self.maccfps = std_tendency.get_maccfps(self)
        # >> if i do it like the above, we can bring the value with just call the attribute
        # >> when we instantiate the class, a certain value would be produced
        self.chem_sim = std_tendency.get_sim(self)
        
        self.datasize = len(self.df) #shape[1] is columns, and shape[0] is rows   
        
        self.df_stdval = std_tendency.get_df_stdval(self)
        self.df_stdsim = std_tendency.get_df_stdsim(self)
        self.df_std_sim_label = std_tendency.get_df_stdval_sim_label(self)
        
    '''
    class BB: 
    
    def __init__(self,name):
        self.name = name 
        self.all = []
        
    def add_some(self,some):  # >> some is from global value, when we want to use this function
        self.all.append(some) 
        # >> in the class method function, we can use class attribute
    '''

    def get_maccfps(self):
        df = self.df
        df['Standard Value'].dropna(axis=0)
        smi = df['Canonical Smiles']
        sd = [Chem.MolFromSmiles(m) for m in smi] 
        maccfps = [MACCSkeys.GenMACCSKeys(m) for m in sd]
        return maccfps
        
    def get_sim(self):
        result = DataStructs.BulkTanimotoSimilarity(self.maccfps[self.chem_num],self.maccfps[:])
        return result 
    
    def get_df_stdval(self):
        df = self.df
        return df['Standard Value'].dropna(axis=0)
    
    
    def get_df_stdsim(self):
        
        df_sim =pd.DataFrame(self.chem_sim)
        
        df_result = pd.concat([self.df_stdval,df_sim],axis=1,ignore_index=True)
        df_result.columns = ['Standard_Value','Similarity']
        return df_result
    
    def get_df_stdval_sim_label(self):
        
        df_result = self.df_stdsim
                
        cluster_label = []

        for row in df_result['Similarity']:
            if row >= 0.8 : 
                cluster_label.append('c4')
            elif row >= 0.6 : 
                cluster_label.append('c3')
            elif row >= 0.4 : 
                cluster_label.append('c2')
            elif row >= 0.2 : 
                cluster_label.append('c1')
            else : 
                cluster_label.append('c0')
        
        df_result['cluster_label'] = cluster_label
        
        for i in range(len(df_result)):
            if df_result.loc[i]['cluster_label'] == 'c4':
                self.c4_stdval.append(df_result.loc[i]['Standard_Value'])
                self.c4_stdval_x.append(i)
    

    def classify_data_as_sim(self):
        
        # sim score between 0.0 - 1.0 
        # 0-0.2  >>> class0
        # 0.2 - 0.4 >>> class1 
        # 0.4 - 0.6 >>> class2 
        # 0.6 - 0.8 >>> class3 
        # 0.8 - 1.0 >>> class4 (interesting label)

        for i in range(len(self.chem_sim)):
            if 0.0 <= self.chem_sim[i] <= 0.2: 
                self.c0.append(self.chem_sim[i])
                self.c0_x.append(i)
            elif 0.2 < self.chem_sim[i] <= 0.4 :
                self.c1.append(self.chem_sim[i])
                self.c1_x.append(i)
            elif 0.4 < self.chem_sim[i] <= 0.6 :
                self.c2.append(self.chem_sim[i])
                self.c2_x.append(i)
            elif 0.6 < self.chem_sim[i] <= 0.8 :
                self.c3.append(self.chem_sim[i])
                self.c3_x.append(i)
            elif 0.8 < self.chem_sim[i] <= 1.0 : 
                self.c4.append(self.chem_sim[i])
                self.c4_x.append(i)
        
    
    def std_scatter_of_total_sim(self):
        
        self.classify_data_as_sim()
    
        plt.scatter(self.c0_x,self.c0,c='blue',label='class0')
        plt.scatter(self.c1_x,self.c1,c='black',label='class1')
        plt.scatter(self.c2_x,self.c2,c='yellow',label='class2')
        plt.scatter(self.c3_x,self.c3,c='green',label='class3')
        plt.scatter(self.c4_x,self.c4,c='red',label='class4')
        plt.show()
        
    
    def std_scatter_of_high_sim(self):
                
        plt.scatter(self.c4_stdval_x,self.c4_stdval,c='red',label='class4')
        

    def c4_variance(self):
        if not self.c4_stdval : 
            return print("--- no value ---\n Please get scatter plot first.")
        else : 
            return np.var(self.c4_stdval)