import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# matplotlib inline
import hddm
import os
from patsy import dmatrix


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
congPosteriorDataCPPN2cTotal = hddm.utils.post_pred_gen(m_congCPPN2cTotal, samples=500, groupby=["subj_idx"])
#congPosteriorData['response']=1
#congPosteriorData.loc[congPosteriorData['rt']<0,'response']=0
#gss = congPosteriorData['rt'].index.to_list()
#gsd = pd.DataFrame(gss, columns=['wftps','notsure','sample'])
#gsd['subj_idx'] = gsd['wftps'].str.replace('wfpt.','',regex=False)
#gst = pd.DataFrame(congPosteriorData.values, columns=['rt'])
#dataPosterior = pd.DataFrame({'subj_idx': gsd['subj_idx'], 'rt': gst['rt']})

################################################################################################
statsPosterior = hddm.utils.post_pred_stats(dataRt,congPosteriorDataCPPN2cTotal)
congPosteriorDataCPPN2cTotal.to_csv('Data/Posteriors/PosteriorCongCPPN2cTotal.csv')
statsPosterior.to_csv('Data/Posteriors/PosteriorCongCPPN2cStatsTotal.csv')