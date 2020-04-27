import pandas as pd 

data_dir = './tutorial/datasets/tox21/'

NR_AhR = 'NR-AhR_wholetraining.csv'
NR_AR = 'NR-AR_wholetraining.csv'
NR_AR_LBD = 'NR-AR-LBD_wholetraining.csv'
NR_Aromatase = 'NR-Aromatase_wholetraining.csv'
NR_ER = 'NR-ER_wholetraining.csv'
NR_ER_LBD = 'NR-ER_wholetraining.csv'
NR_PPAR_gamma = 'NR-PPAR-gamma_wholetraining.csv'
SR_ARE = 'SR-ARE_wholetraining.csv'
SR_ATAD5 = 'SR-ATAD5_wholetraining.csv'
SR_HSE = 'SR-HSE_wholetraining.csv'
SR_MMP = 'SR-MMP_wholetraining.csv'
SR_p53 = 'SR-p53_wholetraining.csv'


df_NR_AhR = pd.read_csv(data_dir+NR_AhR)
df_NR_AR = pd.read_csv(data_dir+NR_AR)
df_NR_AR_LBD = pd.read_csv(data_dir+NR_AR_LBD)
df_NR_Aromatase = pd.read_csv(data_dir+NR_Aromatase)
df_NR_ER = pd.read_csv(data_dir+NR_ER)
df_NR_ER_LBD = pd.read_csv(data_dir+NR_ER_LBD)
df_NR_PPAR_gamma = pd.read_csv(data_dir+NR_PPAR_gamma)
df_SR_ARE = pd.read_csv(data_dir+SR_ARE)
df_SR_ATAD5 = pd.read_csv(data_dir+SR_ATAD5)
df_SR_HSE = pd.read_csv(data_dir+SR_HSE)
df_SR_MMP = pd.read_csv(data_dir+SR_MMP)
df_SR_p53 = pd.read_csv(data_dir+SR_p53)

# total dataframe 
df = pd.concat([df_NR_AhR,df_NR_AR,df_NR_AR_LBD,df_NR_Aromatase,df_NR_ER,df_NR_ER_LBD,
                df_NR_PPAR_gamma,df_SR_ARE,df_SR_ATAD5,df_SR_HSE,df_SR_MMP,df_SR_p53],ignore_index=True)

# check number of duplicated molecules
df.duplicated(subset=['SMILES']).value_counts()

# remove duplicated molecules
df.drop_duplicates(subset=['SMILES'])

# remain only toxic data 
df_toxic = df[df.Label!=0.0]

# remain only toxic SMILES
df_toxic_smiles = df_toxic['SMILES']

# save as csv and smi
df_toxic_smiles.to_csv(data_dir + 'toxic_smiles_wholetraining.csv')
df_toxic_smiles.to_csv(data_dir + 'toxic_smiles_wholetraining.smi')