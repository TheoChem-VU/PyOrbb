from pyorb.orbitals import sfo, mo, adf, dftb  # noqa
from scm import plams
import numpy as np
import matplotlib.pyplot as plt
# from yutility import plot, ensure_list
from TCutility import ensure_list

if __name__ == '__main__':
    rkffile = '../test/orbitals/rkf/BH3NH3.rkf'
    orbs = pyfmo.orbitals.Orbitals(rkffile)
    orbs.rename_fragments(['Donor', 'Acceptor'], ['BH3', 'NH3'])
    to_excel(orbs.sfos['BH3'], orbs.sfos['NH3'])



