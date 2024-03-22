from pyorb.orbitals import sfo, mo, adf, dftb, mixing # noqa
from scm import plams
import numpy as np
import matplotlib.pyplot as plt
# from yutility import plot, ensure_list
from TCutility import ensure_list


if __name__ == '__main__':
    rkffile = '../test/orbitals/rkf/BH3NH3.rkf'
    fig = plt.figure(figsize=(12,12))
    orbitals = pyfmo.orbitals.Orbitals(rkffile)
    MOs = orbitals.mos['HOMO-6':'LUMO']
    SFOsF1 = orbitals.sfos['Donor(HOMO-6)':'Donor(LUMO)']
    SFOsF2 = orbitals.sfos['Acceptor(HOMO-3)':'Acceptor(LUMO)']


    plt.show()

