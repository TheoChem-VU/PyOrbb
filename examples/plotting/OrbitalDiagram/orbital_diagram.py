import matplotlib.pyplot as plt
import pyfmo

plt.figure(figsize=(12,12))
orbitals = pyfmo.orbitals.Orbitals('adf.rkf')
MOs = orbitals.mos['HOMO-3':'LUMO+3']
SFOsF1 = orbitals.sfos['Donor(HOMO-2)':'Donor(LUMO+2)']
SFOsF2 = orbitals.sfos['Acceptor(HOMO-2)':'Acceptor(LUMO+2)']
fragments = list(orbitals.fragments)
pyfmo.plotting.diagram(MOs, orbitals, SFOsF1, SFOsF2)
plt.show()

