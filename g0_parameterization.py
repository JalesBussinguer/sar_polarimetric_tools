# Code to derive the G0 parameters given a set of samples from SAR images

import pandas as pd
import numpy as np

samples = pd.read_csv('D:/thesis_data/VEG_INDICES/samples/stratified/campestre/20m/FC_20170112_20m_patches.csv', sep=',')

dprvi_samples = samples['dprvi_patch_0'].values

def g0_estimator(samples, L):

    m1 = np.mean(samples)
    m2 = np.mean(samples**2)

    m212 = m2/(m1**2)

    a = -2-(L+1)/(L*m212)
    g  = m1*(2+(L+1)/(L*m212))

    return round(a, 6), round(g, 6)

alpha, gamma = g0_estimator(dprvi_samples,  L=4)

