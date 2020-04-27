import argparse
import pandas as pd 

parser = argparse.ArgumentParser(description='argparse')
parser.add_argument('--input','-i',default='./tutorial/datasets/tox21',help='Input smiles dataset')
parser.add_argument('--proteinTarget','-p',required=True, help='target data')
args = parser.parse_args()


def GetSmiles(proteinName):
    listX, listY = [], []
    afile = args.input + '/' + proteinName + '_wholetraining.smiles'
    f = open(afile, "r")
    lines = f.readlines()
    for line in lines:
        splitted = line.split(" ")
        listX.append(splitted[0])#length can vary
        listY.append(float(splitted[1]))
    f.close()
    return listX, listY

smi, label = GetSmiles(args.proteinTarget)

######### option 1 - save as smi and bio file ######### 
"""
df_smi = pd.DataFrame(smi)
df_label = pd.DataFrame(label)
df_smi.to_csv(args.input+'/inputs/'+args.proteinTarget + '_wholetraining.smi', index=False,header=False)
df_label.to_csv(args.input+'/inputs/'+args.proteinTarget + '_wholetraining.dat', index=False,header=False)
"""
######### option 2 - save as csv file with columns ['SMILES','LABEL'] ######### 
df = pd.DataFrame(columns=['SMILES','Label'])
df['SMILES'] = smi 
df['Label'] = label
df.to_csv(args.input+'/'+args.proteinTarget + '_wholetraining.csv',index=False)
# remove duplicated and non-toxic data