import pandas as pd
import matplotlib.pyplot as plt
import numpy as npy
# matplotlib inline
import hddm
import os
from patsy import dmatrix
from arviz import hdi

os.chdir('C:/SHwork/Jneuro/hddmProject/')
os.getcwd()

# N2cCPP mediaton
m_CPPN2c = hddm.load('Data/EEGModelRegressorAcrossCongCPPN2cInteractTotalZscore2')

a_CongCg = m_CPPN2c.nodes_db.node["a_C(cong, Treatment('Cg'))[T.Ig]"]
v_CongCg = m_CPPN2c.nodes_db.node["v_C(cong, Treatment('Cg'))[T.Ig]"]
t_CongCg = m_CPPN2c.nodes_db.node["t_C(cong, Treatment('Cg'))[T.Ig]"]

hddm.analyze.plot_posterior_nodes([a_CongCg], bins=10)
hddm.analyze.plot_posterior_nodes([v_CongCg], bins=10)
hddm.analyze.plot_posterior_nodes([t_CongCg], bins=10)

(a_CongCg.trace()).mean(), (a_CongCg.trace()).std(), (a_CongCg.trace() > 0).mean(), hdi(a_CongCg.trace(), hdi_prob=0.95)
(v_CongCg.trace()).mean(), (v_CongCg.trace()).std(), (v_CongCg.trace() > 0).mean(), hdi(v_CongCg.trace(), hdi_prob=0.95)
(t_CongCg.trace()).mean(), (t_CongCg.trace()).std(), (t_CongCg.trace() > 0).mean(), hdi(t_CongCg.trace(), hdi_prob=0.95)

a_CPPCg = m_CPPN2c.nodes_db.node["a_CPPrslopez"]
v_CPPCg = m_CPPN2c.nodes_db.node["v_CPPrslopez"]
t_CPPCg = m_CPPN2c.nodes_db.node["t_CPPrslopez"]

hddm.analyze.plot_posterior_nodes([a_CPPCg], bins=10)
hddm.analyze.plot_posterior_nodes([v_CPPCg], bins=10)
hddm.analyze.plot_posterior_nodes([t_CPPCg], bins=10)

(a_CPPCg.trace()).mean(), (a_CPPCg.trace()).std(), (a_CPPCg.trace() > 0).mean(), hdi(a_CPPCg.trace(), hdi_prob=0.95)
(v_CPPCg.trace()).mean(), (v_CPPCg.trace()).std(), (v_CPPCg.trace() > 0).mean(), hdi(v_CPPCg.trace(), hdi_prob=0.95)
(t_CPPCg.trace()).mean(), (t_CPPCg.trace()).std(), (t_CPPCg.trace() > 0).mean(), hdi(t_CPPCg.trace(), hdi_prob=0.95)

a_N2cCg = m_CPPN2c.nodes_db.node["a_N2cpeakz"]
v_N2cCg = m_CPPN2c.nodes_db.node["v_N2cpeakz"]
t_N2cCg = m_CPPN2c.nodes_db.node["t_N2cpeakz"]

hddm.analyze.plot_posterior_nodes([a_N2cCg], bins=10)
hddm.analyze.plot_posterior_nodes([v_N2cCg], bins=10)
hddm.analyze.plot_posterior_nodes([t_N2cCg], bins=10)

(a_N2cCg.trace()).mean(), (a_N2cCg.trace()).std(), (a_N2cCg.trace() > 0).mean(), hdi(a_N2cCg.trace(), hdi_prob=0.95)
(v_N2cCg.trace()).mean(), (v_N2cCg.trace()).std(), (v_N2cCg.trace() > 0).mean(), hdi(v_N2cCg.trace(), hdi_prob=0.95)
(t_N2cCg.trace()).mean(), (t_N2cCg.trace()).std(), (t_N2cCg.trace() > 0).mean(), hdi(t_N2cCg.trace(), hdi_prob=0.95)


a_N2cCgInt = m_CPPN2c.nodes_db.node["a_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"]
v_N2cCgInt = m_CPPN2c.nodes_db.node["v_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"]
t_N2cCgInt = m_CPPN2c.nodes_db.node["t_N2cpeakz:C(cong, Treatment('Cg'))[T.Ig]"]

hddm.analyze.plot_posterior_nodes([a_N2cCgInt], bins=10)
hddm.analyze.plot_posterior_nodes([v_N2cCgInt], bins=10)
hddm.analyze.plot_posterior_nodes([t_N2cCgInt], bins=10)

(a_N2cCgInt.trace()).mean(), (a_N2cCgInt.trace()).std(), (a_N2cCgInt.trace() > 0).mean(), hdi(a_N2cCgInt.trace(), hdi_prob=0.95)
(v_N2cCgInt.trace()).mean(), (v_N2cCgInt.trace()).std(), (v_N2cCgInt.trace() > 0).mean(), hdi(v_N2cCgInt.trace(), hdi_prob=0.95)
(t_N2cCgInt.trace()).mean(), (t_N2cCgInt.trace()).std(), (t_N2cCgInt.trace() > 0).mean(), hdi(t_N2cCgInt.trace(), hdi_prob=0.95)

a_CPPslopeCgInt = m_CPPN2c.nodes_db.node["a_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"]
v_CPPslopeCgInt = m_CPPN2c.nodes_db.node["v_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"]
t_CPPslopeCgInt = m_CPPN2c.nodes_db.node["t_CPPrslopez:C(cong, Treatment('Cg'))[T.Ig]"]

hddm.analyze.plot_posterior_nodes([a_CPPslopeCgInt], bins=10)
hddm.analyze.plot_posterior_nodes([v_CPPslopeCgInt], bins=10)
hddm.analyze.plot_posterior_nodes([t_CPPslopeCgInt], bins=10)

(a_CPPslopeCgInt.trace()).mean(), (a_CPPslopeCgInt.trace()).std(), (a_CPPslopeCgInt.trace() > 0).mean(), hdi(a_CPPslopeCgInt.trace(), hdi_prob=0.95)
(v_CPPslopeCgInt.trace()).mean(), (v_CPPslopeCgInt.trace()).std(), (v_CPPslopeCgInt.trace() > 0).mean(), hdi(v_CPPslopeCgInt.trace(), hdi_prob=0.95)
(t_CPPslopeCgInt.trace()).mean(), (t_CPPslopeCgInt.trace()).std(), (t_CPPslopeCgInt.trace() > 0).mean(), hdi(t_CPPslopeCgInt.trace(), hdi_prob=0.95)
#save stats
tvaDF= m_CPPN2c.gen_stats()
tvaDF["bin"]=m_CPPN2c.nodes_db.bin
tvaDF["cong"]=m_CPPN2c.nodes_db.cong
tvaDF["stim"]=m_CPPN2c.nodes_db.stim
tvaDF["type"]=tvaDF.index.str[0]
tvaDF=tvaDF.reset_index()
tvaDF = tvaDF.rename(columns={"index": "subj_idx"})
tvaDF.subj_idx=tvaDF.subj_idx.str[-2:]
tvaDF.subj_idx = tvaDF.subj_idx.str.replace(".", "")

tvaDF.to_csv(r'BehaviourDataAcrossCongN2cCPPInteract.csv')


#save traces
tvaTrac = pd.DataFrame(npy.transpose(npy.array([a_CongCg.trace(), v_CongCg.trace(), t_CongCg.trace(),
                                                a_N2cCg.trace(), v_N2cCg.trace(), t_N2cCg.trace(),
                                                a_CPPCg.trace(), v_CPPCg.trace(), t_CPPCg.trace(),
                                                a_CPPslopeCgInt.trace(), v_CPPslopeCgInt.trace(), t_CPPslopeCgInt.trace(),
                                                a_N2cCgInt.trace(), v_N2cCgInt.trace(), t_N2cCgInt.trace()])),
                   columns=['a_Cong', 'v_Cong', 't_Cong','a_N2c', 'v_N2c', 't_N2c','a_CPPrslope', 'v_CPPrslope','t_CPPrslope',
                            'a_CPPrslopeInt', 'v_CPPrslopeInt','t_CPPrslopeInt','a_N2cInt', 'v_N2cInt','t_N2cInt'])
tvaTrac.to_csv(r'vatTraceAcrossCongN2cCPPInteract.csv')