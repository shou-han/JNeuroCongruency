import pandas as pd
import matplotlib.pyplot as plt
import numpy as npy
# matplotlib inline
import hddm
from patsy import dmatrix
import os
import numpy as np

os.chdir('C:/SHwork/Jneuro/hddmProject/')
os.getcwd()
from kabuki.analyze import gelman_rubin
import multiprocess

from patsy import dmatrix

# assess convergence:
def run_model(id):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as npy
    # matplotlib inline
    import hddm
    from patsy import dmatrix
    import os
    import numpy as np
    os.chdir('C:/SHwork/Jneuro/hddmProject/')
    os.getcwd()
    ## 1 is motion, 2 is colour
    data = hddm.load_csv('Data/HSSMdataRTCSD.csv')
    data = hddm.utils.flip_errors(data)
    data = data[data.run > 2]  # remove counts
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
    data = data[~np.isnan(data.rt)]  # remove nans
    datamean = data.groupby(['subj_idx', 'stim']).mean()
    datamean = datamean.reset_index()

    rtColour = data.rt[data.stim == 2].mean()
    rtMotion = data.rt[data.stim == 1].mean()

    data = data[data.stim == 1]
    data.stim[data.stim == 1] = 'Motion'



    data.cong[data.cong == 2] = "Ig"
    data.cong[data.cong == 1] = "Cg"


    data.rt = data.rt.astype('double')

    dmatrix("C(cong, Treatment('Cg'))", data)

    m_regressor_N2c = hddm.HDDMRegressor(data, ["a ~ 1 + C(cong, Treatment('Cg')) + CPPrslopez + N2cpeakz + CPPrslopez:C(cong, Treatment('Cg')) + N2cpeakz:C(cong, Treatment('Cg'))",
                                                       "v ~ 1 +  C(cong, Treatment('Cg')) + CPPrslopez + N2cpeakz + CPPrslopez:C(cong, Treatment('Cg')) + N2cpeakz:C(cong, Treatment('Cg'))",
                                                       "t ~ 1 +  C(cong, Treatment('Cg')) + CPPrslopez + N2cpeakz + CPPrslopez:C(cong, Treatment('Cg')) + N2cpeakz:C(cong, Treatment('Cg'))"], p_outlier=0.05)
    m_regressor_N2c.find_starting_values()
    m_regressor_N2c.sample(10000, burn=5000,dbname='tracesRegressorAcrossCongCPPN2cInteractTotalZscore%i.db'%id, db='pickle')
    m_regressor_N2c.save('Data/EEGModelRegressorAcrossCongCPPN2cInteractTotalZscore%i'%id)
    return m_regressor_N2c

models=[]
for i in range(5):
    process = multiprocess.Process(target= run_model, args=(i,))
    models.append(process)

for j in models:
    j.start()

for j in models:
    j.join()

models =[]
for i in range(5):
    m = hddm.load('Data/EEGModelRegressorAcrossCongCPPN2cInteractTotalZscore%i'%i)
    models.append(m)
rhats = gelman_rubin(models)
s = pd.DataFrame([rhats])
s = s.transpose()