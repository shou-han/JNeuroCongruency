import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# matplotlib inline
import hddm
import os
from patsy import dmatrix
import seaborn as sns


os.chdir('C:/SHwork/Jneuro/hddmProject/')
os.getcwd()
################################################################################################

data = hddm.load_csv('Data/HSSMdataRTCSD.csv')
data = hddm.utils.flip_errors(data)
data = data[data.run > 2]  # remove counts
data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
data = data[~np.isnan(data.rt)]  # remove nans
data = data[data.stim == 1]
data.stim[data.stim == 1] = 'Motion'

data.cong[data.cong == 2] = "Ig"
data.cong[data.cong == 1] = "Cg"

data.rt = data.rt.astype('double')
dataRt = data[['subj_idx','rt','response']]

################################################################################################
#comparisons
m_congCPPN2cTotal = hddm.load('Data/EEGModelRegressorAcrossCongCPPN2cInteractTotalZscore3')
m_congCPPN2cSim = pd.read_csv('Data/Posteriors/PosteriorCongCPPN2cTotal.csv')
m_congCPPN2cSim['node']=m_congCPPN2cSim['node'].str.replace('^wfpt\.','',regex=True)

subjects = m_congCPPN2cSim['node'].unique()
m_congCPPN2cTotal.data['subj_idx'] =m_congCPPN2cTotal.data['subj_idx'] .astype(str)

num_subjects = len(subjects)
fig, axes = plt.subplots(nrows=11, ncols=4, figsize=(10,6 * num_subjects))

if num_subjects == 1:
    axes = [axes]

axes = axes.flatten()

for i, subject in enumerate(subjects):
    print(f'Subject: {i}')
    for j in range(0, 500, 10):
        #print(f'Subject: {j}')
        subset1 = m_congCPPN2cSim[(m_congCPPN2cSim['node']==subject) & (m_congCPPN2cSim['sample']==j)]
        subset2 = m_congCPPN2cTotal.data[m_congCPPN2cTotal.data['subj_idx'] == subject]
        sns.kdeplot(data=subset1.rt, fill=False, ax=axes[i], label='simulation', color='blue',alpha=0.1)
    sns.kdeplot(data=subset2.rt, fill=False, ax=axes[i], label='data', color='red')
    axes[i].grid(True)
plt.tight_layout()
plt.show()

plt.savefig(f'figure_samples2.png')
