# Code to derive the G0 parameters given a set of samples from SAR images

import numpy as np
import pandas as pd
from scipy.special import gamma
from scipy.optimize import minimize

samples = pd.read_csv('D:/thesis_data/VEG_INDICES/samples/stratified/campestre/20m/FC_20170112_20m_patches.csv', sep=',')

dprvi_samples = samples['dprvi_patch_0'].values

def g0_estimator_moments(samples, L):

    # This function derives the parameters alpha and gamma from a G0 distribution by the moments method
    # This method assumes the Number of Looks (L) as known (for Sentinel-1: L = 4)

    m1 = np.mean(samples)
    m2 = np.mean(samples**2)

    m212 = m2/(m1**2)

    a = -2-(L+1)/(L*m212)
    g  = m1*(2+(L+1)/(L*m212))

    return round(a, 6), round(g, 6)

alpha0, gamma0 = g0_estimator_moments(dprvi_samples,  L=4)

def loglikelihood_estimator(chute_inicial):

    # This function derives the parameters alpha and gamma from a G0 distribution by the maximum likelihhod method
    # The input for this function is the alpha and gamma estimated previously by the moments method
    # This method assumes the Number of Looks (L) as known (for Sentinel-1: L = 4)
    # Samples is an array containing the pixels values of a given region in the SAR image

    p_alpha = chute_inicial[0]
    p_gamma = chute_inicial[1]

    samples = dprvi_samples
    n = len(samples) # Number of samples

    return n * (gamma(4 - p_alpha) - p_alpha * np.log(p_gamma) - gamma(-p_alpha)) + (p_alpha-4) * sum(np.log(p_gamma + samples * 4))

chute_inicial = [alpha0, gamma0]

results = minimize(loglikelihood_estimator, x0=chute_inicial, method='BFGS')
print(f'a = {results.x[0]}, g = {results.x[1]}')