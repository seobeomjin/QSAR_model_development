from rdkit.Chem import MACCSkeys        #MACCS Keys
from rdkit.Chem.AtomPairs import Pairs, Torsions  #Atom Pairs and Torsions      
from rdkit import Chem, DataStructs

import numpy as np
import pandas as pd
import os 


class Descriptor() : 

    def __init__(self, smi):
        self.smi = smi
        self.sd = [Chem.MolFromSmiles(m) for m in smi]
        
        #output
        self.maccfps = Descriptor.Calc_maccfps(self)
        self.maccfps_bit = Descriptor.Calc_maccfps_bit(self) 
    
    #MACCS keys 
    def Calc_maccfps(self):
        maccfps = [MACCSkeys.GenMACCSKeys(m) for m in self.sd]
        return maccfps
    
    def Calc_maccfps_bit(self):
        np_fps = []
    
        if self.maccfps == None :
            print("Empty")
        else:
            for fp in self.maccfps :
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp,arr)
                np_fps.append(arr)

        return np_fps

    #Atom Pairs
    def Calc_AtomPairs(self): #type error
        pairFps = [Pairs.GetAtomPairFingerprint(x) for x in self.sd]
        return pairFps
    
    def Calc_AtomPairs_Int(self): #type error
        pairFps_int = [Pairs.GetAtomPairFingerprintAsIntVect(x) for x in self.sd]
        return pairFps_int

    def Calc_AtomPairs_Bit(self):
        pairFps_bit = [Pairs.GetAtomPairFingerprintAsBitVect(x) for x in self.sd]
        return pairFps_bit

    #Torsions 
    def Calc_Torsions(self):
        tts = [Torsions.GetTopologicalTorsionFingerprintAsIntVect(x) for x in self.sd]
        return tts 

    #E3FP

#path for deleted Nan Value on Standard Value
wdir_dltnan = 'C:/jupyter_devel/kist-europe/QSAR/AOP_data/dltnan/'

wdir = 'C:/jupyter_devel/kist-europe/QSAR/AOP_data/dltnan/'
csv = 'final_ic50_dltnan.csv'

df = pd.read_csv(wdir + csv)
smi = df['Canonical Smiles']

descriptor = Descriptor(smi)
Vect = descriptor.Calc_AtomPairs_Bit() 
fp_df = pd.DataFrame(Vect,columns=list(range(len(Vect[0]))))

df_id = df[['Molecule','Standard Value']]
#df_id = df[['Molecule ChEMBL ID','Standard Value']] << for androgen dataset 

#Data table preparation
df_label_fp = pd.concat([df_id, fp_df], axis=1)
df_label_fp.to_csv(wdir+'5AR_check-atompair-stdval.csv',index=False)


