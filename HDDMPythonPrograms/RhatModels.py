import pandas as pd
import matplotlib.pyplot as plt
import numpy as npy
# matplotlib inline
import hddm
import os
from patsy import dmatrix
from arviz import hdi
from kabuki.analyze import gelman_rubin

models =[]
for i in range(5):
    m_N2c = hddm.load('Data/EEGModelRegressorAcrossN2cZscore%i'%i)
    models.append(m_N2c)
rhatsN2c = gelman_rubin(models)
rhatsN2c["a_N2cpeakz"], rhatsN2c["v_N2cpeakz"], rhatsN2c["t_N2cpeakz"]
s_N2c = pd.DataFrame([rhatsN2c])
s_N2c = s_N2c.transpose()

models =[]
for i in range(5):
    m_CPPN2c = hddm.load('Data/EEGModelRegressorAcrossWithinCPPN2cZscore%i'%i)
    models.append(m_CPPN2c)
rhatsCPPN2c = gelman_rubin(models)
rhatsCPPN2c["a_N2cpeakz"], rhatsCPPN2c["v_N2cpeakz"], rhatsCPPN2c["t_N2cpeakz"],rhatsCPPN2c["a_CPPrslopez"], rhatsCPPN2c["v_CPPrslopez"], rhatsCPPN2c["t_CPPrslopez"]
s_CPPN2c = pd.DataFrame([rhatsCPPN2c])
s_CPPN2c = s_CPPN2c.transpose()

models =[]
for i in range(5):
    m_congN2c = hddm.load('Data/EEGModelRegressorAcrosscongN2cZscore%i'%i)
    models.append(m_congN2c)
rhatscongN2c = gelman_rubin(models)
rhatscongN2c["a_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongN2c["v_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongN2c["t_C(cong, Treatment('Cg'))[T.Ig]"],rhatscongN2c["a_N2cpeakz"], rhatscongN2c["v_N2cpeakz"], rhatscongN2c["t_N2cpeakz"]
s_congN2c = pd.DataFrame([rhatscongN2c])
s_congN2c = s_congN2c.transpose()

models =[]
for i in range(5):
    m_cong = hddm.load('Data/EEGModelRegressorAcrosscongWithinZscore%i'%i)
    models.append(m_cong)
rhatscong = gelman_rubin(models)
rhatscong["a_C(cong, Treatment('Cg'))[T.Ig]"], rhatscong["v_C(cong, Treatment('Cg'))[T.Ig]"], rhatscong["t_C(cong, Treatment('Cg'))[T.Ig]"]
s_cong = pd.DataFrame([rhatscong])
s_cong = s_cong.transpose()


models =[]
for i in range(5):
    m_congCPP = hddm.load('Data/EEGModelRegressorAcrosscongCPPZscore%i'%i)
    models.append(m_congCPP)
rhatscongCPP = gelman_rubin(models)
rhatscongCPP["a_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPP["v_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPP["t_C(cong, Treatment('Cg'))[T.Ig]"],rhatscongCPP["a_CPPrslopez"], rhatscongCPP["v_CPPrslopez"], rhatscongCPP["t_CPPrslopez"]

s_congCPP = pd.DataFrame([rhatscongCPP])
s_congCPP = s_congCPP.transpose()



models =[]
for i in range(5):
    m_congCPPN2c = hddm.load('Data/EEGModelRegressorAcrossCongCPPN2cZscore%i'%i)
    models.append(m_congCPPN2c)
rhatscongCPPN2c = gelman_rubin(models)
s_congCPPN2c= pd.DataFrame([rhatscongCPPN2c])
s_congCPPN2c = s_congCPPN2c.transpose()
(rhatscongCPPN2c["a_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2c["v_C(cong, Treatment('Cg'))[T.Ig]"],  rhatscongCPPN2c["t_C(cong, Treatment('Cg'))[T.Ig]"],
 rhatscongCPPN2c["a_N2cpeakz"], rhatscongCPPN2c["v_N2cpeakz"], rhatscongCPPN2c["t_N2cpeakz"],
 rhatscongCPPN2c["a_CPPrslopez"], rhatscongCPPN2c["v_CPPrslopez"], rhatscongCPPN2c["t_CPPrslopez"]

 )



models =[]
for i in range(5):
    m_congCPPN2cTotal = hddm.load('Data/EEGModelRegressorAcrossCongCPPN2cInteractTotalZscore%i'%i)
    models.append(m_congCPPN2cTotal)
rhatscongCPPN2cTotal = gelman_rubin(models)
s_congCPPN2cTotal= pd.DataFrame([rhatscongCPPN2cTotal])
s_congCPPN2cTotal = s_congCPPN2cTotal.transpose()
(rhatscongCPPN2cTotal["a_C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2cTotal["v_C(cong, Treatment('Cg'))[T.Ig]"],  rhatscongCPPN2cTotal["t_C(cong, Treatment('Cg'))[T.Ig]"],
 rhatscongCPPN2cTotal["a_N2cpeakz"], rhatscongCPPN2cTotal["v_N2cpeakz"], rhatscongCPPN2cTotal["t_N2cpeakz"],
 rhatscongCPPN2cTotal["a_CPPrslopez"], rhatscongCPPN2cTotal["v_CPPrslopez"], rhatscongCPPN2cTotal["t_CPPrslopez"],
 rhatscongCPPN2cTotal["a_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2cTotal["v_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2cTotal["t_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"],
 rhatscongCPPN2cTotal["a_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2cTotal["v_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"], rhatscongCPPN2cTotal["t_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"]
 )