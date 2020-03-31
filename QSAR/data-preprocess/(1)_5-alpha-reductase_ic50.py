# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 00:58:46 2019

@author: phdgil
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#wdir is directory where your data is located.
wdir ='/home/beomjin/kist-europe/QSAR_model_development/AOP_data/'
xlsx = '5-alpha-reductase1_finalIC50.xlsx' #'ic50.xlsx'

df = pd.read_excel(wdir+xlsx)

#Remove if data is 'outside typical range'
df = df[df['Data Validity Comment'].isna()]

#Remove certain assay
remove_assay = ['CHEMBL810002', 'CHEMBL810006']
for r in remove_assay:
    df = df[df['Assay']!=r]

#Checking duplicated CHEMBL ID for compounds
df['chemical_duplicated'] = 0
molist = df['Molecule'].tolist()
molist = list(set(molist))
for m in molist:
    if df[df['Molecule']==m].shape[0]!=1: 
        df.loc[df['Molecule']==m,'chemical_duplicated']=1

#labeling active, inactive
df['label']=0
df.loc[df['Standard Value']<100,'label']=1

#Data selection (Same Chembl ID and same label)
df_duple = df[df['chemical_duplicated']==1]

molist = df_duple['Molecule'].tolist()
molist = list(set(molist))
for m in molist:
    if df_duple[df_duple['Molecule']==m]['label'].mean() == 0:
        df = df.drop(df[df['Molecule']==m].index.values[1:])
        df_duple = df_duple.drop(df_duple[df_duple['Molecule']==m].index.values)
    elif df_duple[df_duple['Molecule']==m]['label'].mean() == 1:
        df = df.drop(df[df['Molecule']==m].index.values[1:])
        df_duple = df_duple.drop(df_duple[df_duple['Molecule']==m].index.values)
    else:
        print('Label mismatch: ', m)

df = df.drop(df[df['Molecule']=='CHEMBL1201841'].index.values[1:])
df = df.drop(df[df['Molecule']=='CHEMBL1642918'].index.values)
df = df.drop(['chemical_duplicated'], axis=1)

#Checking duplicated MW
df['mw_duplicated']=0
mwlist = df['Molecular Weight'].tolist()
mwlist = list(set(mwlist))
for m in mwlist:
    if df[df['Molecular Weight']==m].shape[0]!=1: ### if shape[0] of df which the value of MW is same with m is not 1(it means true?) ,,?
        df.loc[df['Molecular Weight']==m, 'mw_duplicated']=1

#Each chemical structures were visually inspected to remove duplicated chemical structures.
#remove_redundant = [177, 229, 90, 359, 77, 117, 263, 165, 234, 348, 361, 119,
#                    276, 473, 19, 103, 275, 95, 216, 371, 399, 401, 20, 41, 63,
#                    179, 21, 9, 280, 118, 169, 530, 80, 62, 199, 164, 197, 443,
#                    363, 59, 116, 404, 128, 241, 282, 434, 449, 127, 594, 360,
#                    156, 447, 189, 313, 267, 113, 167, 545, 230, 214, 382,
#                    10, 97, 279, 474, 91, 441, 383, 292, 40]
### it looks like indices of chemical structures. number is spreading to 545,, but the actual data is numbered to 288

#df = df.drop(remove_redundant)

df.to_csv(wdir+'final_ic50.csv',index=False)
'''
#smiles to sdf (For structure examination)

import openbabel
import os

babel_convert=openbabel.OBConversion()
babel_convert.SetInAndOutFormats('smi','mol')

smidf = df[df['mw_duplicated']==1]
smindex = smidf.index.values

OBmol = openbabel.OBMol()
for i in smindex:
    babel_convert.ReadString(OBmol, smidf.loc[i,'Canonical Smiles'])
    subdir = 'mol/' + str(smidf.loc[i,'Molecular Weight']) +'/'
    if not os.path.exists(wdir+subdir):
        os.makedirs(wdir+subdir)
    babel_convert.WriteFile(OBmol, wdir+subdir+str(i)+'_'+smidf.loc[i,'Molecule']+'.mol') 

#List of structures to be removed (duplicated structures)
#mw: 331.46
df = df.drop(df[df['Molecule']=='CHEMBL300446'].index.values) #177
df = df.drop(df[df['Molecule']=='CHEMBL297524'].index.values) #229

#mw: 332.36
df = df.drop(df[df['Molecule']=='CHEMBL382649'].index.values) #90

#mw: 343.51
df = df.drop(df[df['Molecule']=='CHEMBL311577'].index.values) #359

#mw: 357.54
df = df.drop(df[df['Molecule']=='CHEMBL420020'].index.values) #77

#mw: 358.53
df = df.drop(df[df['Molecule']=='CHEMBL136863'].index.values) #117
df = df.drop(df[df['Molecule']=='CHEMBL135313'].index.values) #263

#mw: 360.5
df = df.drop(df[df['Molecule']=='CHEMBL307626'].index.values) #165

#mw: 370.54
df = df.drop(df[df['Molecule']=='CHEMBL3706583'].index.values) #234
df = df.drop(df[df['Molecule']=='CHEMBL3706586'].index.values) #348

#mw: 371.57
df = df.drop(df[df['Molecule']=='CHEMBL80160'].index.values) #361

#mw: 372.55
df = df.drop(df[df['Molecule']=='CHEMBL76192'].index.values) #119
df = df.drop(df[df['Molecule']=='CHEMBL76648'].index.values) #276

#mw: 384.56
df = df.drop(df[df['Molecule']=='CHEMBL309261'].index.values) #473
df = df.drop(df[df['Molecule']=='CHEMBL3706585'].index.values) #19
df = df.drop(df[df['Molecule']=='CHEMBL306507'].index.values) #103
df = df.drop(df[df['Molecule']=='CHEMBL77129'].index.values) #275

#mw: 385.55
df = df.drop(df[df['Molecule']=='CHEMBL310710'].index.values) #95

#mw: 386.54
df = df.drop(df[df['Molecule']=='CHEMBL24291'].index.values) #216

#mw: 386.58
df = df.drop(df[df['Molecule']=='CHEMBL78275'].index.values) #371
df = df.drop(df[df['Molecule']=='CHEMBL77765'].index.values) #399
df = df.drop(df[df['Molecule']=='CHEMBL76892'].index.values) #401

#mw: 388.55
df = df.drop(df[df['Molecule']=='CHEMBL78064'].index.values) #20

#mw: 392.54
df = df.drop(df[df['Molecule']=='CHEMBL280155'].index.values) #41
df = df.drop(df[df['Molecule']=='CHEMBL324210'].index.values) #63

#mw: 397.56
df = df.drop(df[df['Molecule']=='CHEMBL305932'].index.values) #179

#mw: 398.59
df = df.drop(df[df['Molecule']=='CHEMBL80669'].index.values) #21

#mw: 400.61
df = df.drop(df[df['Molecule']=='CHEMBL75320'].index.values) #9
df = df.drop(df[df['Molecule']=='CHEMBL77144'].index.values) #280
df = df.drop(df[df['Molecule']=='CHEMBL137691'].index.values) #118
df = df.drop(df[df['Molecule']=='CHEMBL308532'].index.values) #169

#mw: 401.68
df = df.drop(df[df['Molecule']=='CHEMBL340027'].index.values) #530

#mw: 402.6
df = df.drop(df[df['Molecule']=='CHEMBL73816'].index.values) #80

#mw: 406.57
df = df.drop(df[df['Molecule']=='CHEMBL111407'].index.values) #62
df = df.drop(df[df['Molecule']=='CHEMBL135780'].index.values) #199

#mw: 407
df = df.drop(df[df['Molecule']=='CHEMBL2282660'].index.values) #164

#mw: 414.63
df = df.drop(df[df['Molecule']=='CHEMBL80163'].index.values) #197
df = df.drop(df[df['Molecule']=='CHEMBL76536'].index.values) #443

#mw: 416.57
df = df.drop(df[df['Molecule']=='CHEMBL24465'].index.values) #363

#mw: 416.65
df = df.drop(df[df['Molecule']=='CHEMBL169621'].index.values) #59

#mw: 428.66
df = df.drop(df[df['Molecule']=='CHEMBL77378'].index.values) #116

#mw: 429.73
df = df.drop(df[df['Molecule']=='CHEMBL340305'].index.values) #404

#mw: 436.43
df = df.drop(df[df['Molecule']=='CHEMBL311045'].index.values) #128

#mw: 442.6
df = df.drop(df[df['Molecule']=='CHEMBL107654'].index.values) #241

#mw: 448.65
df = df.drop(df[df['Molecule']=='CHEMBL88494'].index.values) #282

#mw: 450.67
df = df.drop(df[df['Molecule']=='CHEMBL2282782'].index.values) #434
df = df.drop(df[df['Molecule']=='CHEMBL77024'].index.values) #449

#mw: 451.45
df = df.drop(df[df['Molecule']=='CHEMBL307181'].index.values) #127

#mw: 451.65
df = df.drop(df[df['Molecule']=='CHEMBL2282783'].index.values) #594

#mw: 456.72
df = df.drop(df[df['Molecule']=='CHEMBL421696'].index.values) #360

#mw: 460.54
df = df.drop(df[df['Molecule']=='CHEMBL110288'].index.values) #156

#mw: 462.68
df = df.drop(df[df['Molecule']=='CHEMBL76671'].index.values) #447

#mw: 468.64
df = df.drop(df[df['Molecule']=='CHEMBL322749'].index.values) #189
df = df.drop(df[df['Molecule']=='CHEMBL321442'].index.values) #313

#mw: 477.78
df = df.drop(df[df['Molecule']=='CHEMBL339912'].index.values) #267

#mw: 482.67
df = df.drop(df[df['Molecule']=='CHEMBL323996'].index.values) #113
df = df.drop(df[df['Molecule']=='CHEMBL24088'].index.values) #167

#mw: 483.66
df = df.drop(df[df['Molecule']=='CHEMBL312531'].index.values) #545

#mw: 494.76
df = df.drop(df[df['Molecule']=='CHEMBL77328'].index.values) #230

#mw: 496.7
df = df.drop(df[df['Molecule']=='CHEMBL436215'].index.values) #214
df = df.drop(df[df['Molecule']=='CHEMBL306554'].index.values) #382

#mw: 498.45
df = df.drop(df[df['Molecule']=='CHEMBL76803'].index.values) #10

#mw: 503.77
df = df.drop(df[df['Molecule']=='CHEMBL77124'].index.values) #97

#mw: 518.65
df = df.drop(df[df['Molecule']=='CHEMBL306722'].index.values) #279

#mw: 522.73
df = df.drop(df[df['Molecule']=='CHEMBL309585'].index.values) #474

#mw: 528.54
df = df.drop(df[df['Molecule']=='CHEMBL89240'].index.values) #91
df = df.drop(df[df['Molecule']=='CHEMBL108559'].index.values) #441

#mw: 551.56
df = df.drop(df[df['Molecule']=='CHEMBL412425'].index.values) #383

#mw: 558.77
df = df.drop(df[df['Molecule']=='CHEMBL308258'].index.values) #292

#mw: 592.75
df = df.drop(df[df['Molecule']=='CHEMBL314069'].index.values) #40
'''