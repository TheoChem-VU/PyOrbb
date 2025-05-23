import pyfmo
import tcutility
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools as it

p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.results/'
orbs = pyfmo.orbitals2.objects.Orbitals(p + 'adf.rkf')
res = tcutility.results.read(p)
print(res.properties.energy.orbint.total)

f1, f2 = orbs.fragments
sfos1, sfos2 = orbs.sfos.get_fragment_sfos(f1), orbs.sfos.get_fragment_sfos(f2)
sfos1_occ, sfos2_occ = np.array([sfo.occupied for sfo in sfos1]).reshape(1, -1), np.array([sfo.occupied for sfo in sfos2]).reshape(-1, 1)

E1, E2 = np.array([sfo.energy for sfo in sfos1]).reshape(1, -1), np.array([sfo.energy for sfo in sfos2]).reshape(-1, 1)
dE = abs(E2 - E1)
S = orbs.data.matrices.overlap.total
S = orbs.overlap_matrix(sfos1, sfos2)

occ_virt_mask = np.logical_or(sfos1_occ, sfos2_occ)
oi = -S**2 / dE * occ_virt_mask

idxs = np.argsort(oi, axis=None)
cs = oi.flatten()
cs = np.sort(cs)
cs = np.cumsum(cs) / sum(cs)

threshold = 0.7
chosen_interactions = []
for i, x in enumerate(cs):
	if x <= threshold:
		idx = np.unravel_index(idxs[i], oi.shape)
		chosen_interactions.append((*idx[::-1], oi[idx], x))
	else:
		break

def argNmax(arr, index):
	return np.argsort(arr)[-index-1]

assoc_mos = {}
G_temp = nx.Graph()
for idx1, idx2, strength, cumsum_strength in chosen_interactions:
	sfo1, sfo2 = sfos1[idx1], sfos2[idx2]
	mo_occ_mask = np.array([mo.occupied for mo in orbs.mos])
	mo_virt_mask = 1 - mo_occ_mask
	N_occ_sfos = sfo1.occupied + sfo2.occupied

	contr1 = abs(orbs.data.matrices.mulliken_contribution.total[:, sfo1.index-1])
	contr2 = abs(orbs.data.matrices.mulliken_contribution.total[:, sfo2.index-1])
	contr = contr1 * contr2

	if N_occ_sfos == 0:
		mo1 = orbs.mos.orbitals[argNmax(contr * mo_occ_mask, 0)]
		mo2 = orbs.mos.orbitals[argNmax(contr * mo_occ_mask, 1)]
	if N_occ_sfos == 1:
		mo1 = orbs.mos.orbitals[argNmax(contr * mo_occ_mask, 0)]
		mo2 = orbs.mos.orbitals[argNmax(contr * mo_virt_mask, 0)]
	if N_occ_sfos == 2:
		mo1 = orbs.mos.orbitals[argNmax(contr * mo_virt_mask, 0)]
		mo2 = orbs.mos.orbitals[argNmax(contr * mo_virt_mask, 1)]

	assoc_mos[tuple(sorted([sfo1, sfo2], key=lambda sfo: sfo.energy))] = (mo1, mo2)

	G_temp.add_node(sfo1, )
	G_temp.add_node(sfo2, pos=(1, sfo2.energy))
	# G_temp.add_node(mo1, pos=(0, mo1.energy))
	# G_temp.add_node(mo2, pos=(0, mo2.energy))
	
	G_temp.add_edge(sfo1, sfo2, weight=strength)
	# G_temp.add_edge(sfo1, mo2, weight=0)
	# G_temp.add_edge(sfo2, mo1, weight=0)
	# G_temp.add_edge(sfo2, mo2, weight=0)


components = [G_temp.subgraph(H) for H in nx.connected_components(G_temp)]
components = sorted(components, key=lambda H: sum(H.get_edge_data(*e, default={'weight': 0})['weight'] for e in H.edges()))

for H in components:
	G_final = nx.Graph()
	sfos = H.nodes()
	sfos1 = [sfo for sfo in sfos if sfo.fragment == f1]
	sfos2 = [sfo for sfo in sfos if sfo.fragment == f2]

	plt.figure()
	# pos = nx.spring_layout(H)

	for (sfo1, sfo2) in it.product(sfos1, sfos2):
		sorted_ = tuple(sorted([sfo1, sfo2], key=lambda sfo: sfo.energy))
		if sorted_ not in assoc_mos:
			continue
			# print(sorted_[0], *assoc_mos[sorted_], sorted_[1])
		mo1, mo2 = assoc_mos[sorted_]
		G_final.add_node(sfo1, pos=(-1, sfo1.energy))
		G_final.add_node(sfo2, pos=( 1, sfo2.energy))
		G_final.add_node(mo1,  pos=( 0,  mo1.energy))
		G_final.add_node(mo2,  pos=( 0,  mo2.energy))

		nx.draw_networkx_nodes(G_final, G_final.nodes(data='pos'))
		nx.draw_networkx_labels(G_final, G_final.nodes(data='pos'))

		G_final.add_edge(sfo1, mo1)
		G_final.add_edge(sfo1, mo2)
		G_final.add_edge(sfo2, mo1)
		G_final.add_edge(sfo2, mo2)
		nx.draw_networkx_edges(G_final, G_final.nodes(data='pos'), edgelist=[(sfo1, mo1)], alpha=np.clip(abs(sfo1.mulliken_contribution(mo1)), 0, 1))
		nx.draw_networkx_edges(G_final, G_final.nodes(data='pos'), edgelist=[(sfo1, mo2)], alpha=np.clip(abs(sfo1.mulliken_contribution(mo2)), 0, 1))
		nx.draw_networkx_edges(G_final, G_final.nodes(data='pos'), edgelist=[(sfo2, mo1)], alpha=np.clip(abs(sfo2.mulliken_contribution(mo1)), 0, 1))
		nx.draw_networkx_edges(G_final, G_final.nodes(data='pos'), edgelist=[(sfo2, mo2)], alpha=np.clip(abs(sfo2.mulliken_contribution(mo2)), 0, 1))


	weights = nx.get_edge_attributes(H, 'weight')
	plt.title(sum(weights.values()))
	weights = {k: round(v, 3) for k, v in weights.items()}

plt.figure()
plt.plot(cs)
plt.hlines(threshold, 0, len(cs), colors='r', linestyle='dashed')
plt.show()
