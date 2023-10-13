import numpy as np
from pyorb.orbitals import sfo


def get_alpha(sfos, mos):
    '''
    Calculate alpha, which is the sum of squared coefficients of each sfo in each mo
    '''

    alphas = []
    for mo in mos:
        coeffs = np.array(mo.get_coeff(sfo))/sum(abs(mo.coeffs**2))
        alphas.append(sum(coeffs**2))
    return np.array(alphas)


def predict_interactions(sfos1, mos1, sfos2, mos2, use_mask=True):
    alphas1 = get_alpha(sfos1, mos1).reshape(-1, 1)
    alphas2 = get_alpha(sfos2, mos2).reshape(-1, 1)

    energy1 = np.array([mo.energy for mo in mos1]).reshape(-1, 1)
    energy2 = np.array([mo.energy for mo in mos2]).reshape(-1, 1)

    S = alphas1 * alphas2.T
    dE = abs(energy1 - energy2.T)
    oi = S/dE
    if use_mask:
        mask = sfo.occ_virt_mask(mos1, mos2)
        oi[np.logical_not(mask)] = None

    return oi
