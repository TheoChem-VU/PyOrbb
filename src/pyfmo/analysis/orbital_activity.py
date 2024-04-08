import numpy as np


def mos(orbs):
    '''
    Get activity terms of MOs.
    MO activity is defined as the magnitude of the difference in the Mulliken population and the occupation.
    A bigger difference should indicate that the MO originates from stronger charge-transfer interaction.
    '''
    P = orbs.mulliken_population_matrix()
    Pmarginal = np.sum(abs(P), axis=0)
    activity = abs(Pmarginal - np.array([mo.occupation for mo in orbs.mos]))
    return activity


def sfos(orbs):
    '''
    Get activity terms of SFOs.
    SFO activity is defined as the sum of the absolute Mulliken-contributions minus 1.
    The Mulliken-contributions for an SFO should sum up to 1. 
    However, since Mulliken-contributions can be positive or negative, 
    it is more interesting to look at the sum of the absolute values.
    '''
    C = orbs.mulliken_contribution_matrix()
    Cabs = abs(C)
    Cabs_marginal = np.sum(Cabs, axis=0)
    return Cabs_marginal - 1
