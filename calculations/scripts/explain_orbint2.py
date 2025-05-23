import pyfmo
import tcutility
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools as it


def scree(arr, thresh):
	# arr = (arr - arr.min()) / (arr.max() - arr.min())
	arr = abs(arr)
	arr = arr / arr.sum()
	arr = np.sort(arr)[::-1]
	cs = np.cumsum(arr)
	# plt.plot(cs)
	# plt.show()
	return sum(cs < thresh)

p = '/Users/yumanhordijk/PhD/Projects/RadicalAdditionASMEDA/data/DFT/TS_X_O/PyFrag_OLYP_TZ2P/Step.00035/complex/'
orbs = pyfmo.orbitals2.objects.Orbitals(p + 'adf.rkf')
res = tcutility.results.read(p)
print(res.properties.energy.orbint.total)

dE = []
dE_dict = {}
ncomponents = []
MO_energies = np.array([mo.energy for mo in orbs.mos])
MO_occ = np.array([mo.occupied for mo in orbs.mos])
for sfo in orbs.sfos:
	contr = abs(orbs.data.matrices.mulliken_contribution.total[:, sfo.index - 1])
	print(sfo, sfo.energy, sum(MO_energies * contr))
	dE.append((MO_energies * contr) / sum(MO_energies * contr))
	dE_dict[sfo] = sum(MO_energies * contr) - sfo.energy
	ncomponents.append(scree(MO_energies * contr, 0.75))

dE = np.array(dE)
print(np.sum(dE, axis=1))
# plt.imshow(dE, origin='lower')
plt.plot(ncomponents)
plt.show()
