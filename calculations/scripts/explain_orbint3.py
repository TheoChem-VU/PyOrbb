import pyfmo
import tcutility
import numpy as np
import matplotlib.pyplot as plt
# import networkx as nx
import itertools as it
# import tcviewer


def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in tcutility.ensure_list(sfos1):
        ret.append([])
        for sfo2 in tcutility.ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


def argNmax(arr, N):
    '''
    Get the Nth maximum element of an array.
    '''
    return np.argsort(arr, axis=None)[-N-1]


def _get_orbint(orbs, n_pairs=20):
    sfos = {}
    sfo_occ = {}
    sfo_ene = {}
    for i, frag in enumerate(orbs.sfos.fragments):
        sfos[frag] = [sfo for sfo in orbs.sfos if sfo.fragment_unique == frag]
        sfo_occ[frag] = np.array([sfo.occupation for sfo in sfos[frag]]).reshape(-1, 1)
        sfo_ene[frag] = np.array([sfo.energy for sfo in sfos[frag]]).reshape(-1, 1)

    mos = list(orbs.mos)
    mos = list(sorted(mos, key=lambda mo: mo.energy))
    occmo = np.array([mo.occupied for mo in mos])

    ret = []
    oi_approx_total = 0
    for i, frag in enumerate(orbs.sfos.fragments):
        for frag2 in orbs.sfos.fragments[i+1:]:
            occ_virt_mask = np.logical_xor(sfo_occ[frag], sfo_occ[frag2].T)
            S = overlap_mat(sfos[frag], sfos[frag2])
            dE = abs(sfo_ene[frag] - sfo_ene[frag2].T)
            nogap_mask = dE != 0
            total_mask = occ_virt_mask * nogap_mask
            dE += (1 - nogap_mask)
            oi = -S**2 / dE * total_mask
            oi_approx_total += oi.sum()
            for i in range(n_pairs):
                best = np.unravel_index(argNmax(-oi, i), oi.shape)
                sfo1, sfo2  = sfos[frag][best[0]], sfos[frag2][best[1]]

                sfo1_contr = np.array([sfo1.mulliken_contribution(mo) for mo in mos])
                sfo2_contr = np.array([sfo2.mulliken_contribution(mo) for mo in mos])

                occ_contrs = abs(sfo1_contr * sfo2_contr) * occmo
                virt_contrs = abs(sfo2_contr * sfo1_contr) * (1-occmo)

                occ_mo = mos[np.argmax(occ_contrs)]
                virt_mo = mos[np.argmax(virt_contrs)]

                ret.append((oi[best], sfo1, sfo2, occ_mo, virt_mo))

    Eoi = orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503

    ret = sorted(ret, key=lambda row: row[0])
    ret = ret[:n_pairs]
    ret = [(row[0] * Eoi / oi_approx_total, row[0] / oi_approx_total, *row[1:]) for row in ret]
    return ret


def _get_pauli(orbs, n_pairs=20):
    sfos = {}
    sfo_occ = {}
    sfo_ene = {}
    for i, frag in enumerate(orbs.sfos.fragments):
        sfos[frag] = [sfo for sfo in orbs.sfos if sfo.fragment_unique == frag]
        sfo_occ[frag] = np.array([sfo.occupation for sfo in sfos[frag]]).reshape(-1, 1)

    mos = [mo for mo in orbs.mos if mo.occupied]

    ret = []
    pauli_approx_total = 0
    for i, frag in enumerate(orbs.sfos.fragments):
        for frag2 in orbs.sfos.fragments[i+1:]:
            occ_occ_mask = np.logical_and(sfo_occ[frag], sfo_occ[frag2].T)
            S = overlap_mat(sfos[frag], sfos[frag2])
            pauli = S**2 * occ_occ_mask
            pauli_approx_total += pauli.sum()
            for i in range(n_pairs):
                best = np.unravel_index(argNmax(pauli, i), pauli.shape)
                sfo1, sfo2  = sfos[frag][best[0]], sfos[frag2][best[1]]

                sfo1_contr = np.array([sfo1.mulliken_contribution(mo) for mo in mos])
                sfo2_contr = np.array([sfo2.mulliken_contribution(mo) for mo in mos])

                occ_contrs = abs(sfo1_contr * sfo2_contr)

                occ_mo1 = mos[argNmax(occ_contrs, 0)]
                occ_mo2 = mos[argNmax(occ_contrs, 1)]

                ret.append((pauli[best], sfo1, sfo2, occ_mo1, occ_mo2))

    E_pauli = orbs.reader.read('Energy', 'Pauli Total') * 627.503

    ret = sorted(ret, key=lambda row: -row[0])
    ret = ret[:n_pairs]
    ret = [(row[0] * E_pauli / pauli_approx_total, row[0] / pauli_approx_total, *row[1:]) for row in ret]

    return ret


# p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.results/'
p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/TransitionState/DielsAlder.results/'
# p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/CoordinationBondFeCO4CO/FeCO4CO.results/'
# p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/CoordinationBondFeCO4CH4/FeCO4CH4.results/'
# p = '/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/HydrogenBond/GuanineCytosine.results/'
orbs = pyfmo.orbitals2.objects.Orbitals(p + 'adf.rkf')
res = tcutility.results.read(p)


pairs = _get_orbint(orbs, 20)
# mos, sfos, connections  = [], [], []
for i, (dE, explanation, sfo1, sfo2, mo1, mo2) in enumerate(pairs):
    print(dE, explanation, sfo1, sfo2, mo1, mo2)
    # mos.append([mo1, mo2])
    # sfos.append([sfos1, sfos2])
    # connections.append([(sfo1, mo1), (sfo1, mo2), (sfo2, mo1), (sfo2, mo2)])

    # with tcviewer.Screen() as scr:
    #     with scr.add_molscene() as scene:
    #         cub1 = sfo1.cube_file()
    #         cub2 = sfo2.cube_file()
    #         S = cub1.copy()
    #         S.values *= cub2.values

    #         print(sfo1 @ sfo2)
    #         # scene.draw_isosurface(cub1, -.03, opacity=.25, color=[1, 0, 0])
    #         # scene.draw_isosurface(cub1,  .03, opacity=.25, color=[0, 0, 1])
    #         # scene.draw_molecule(sfo1.molecule)
    #         # scene.draw_isosurface(cub2, -.03, opacity=.25, color=[0, 1, 1])
    #         # scene.draw_isosurface(cub2,  .03, opacity=.25, color=[1, .5, 0])
    #         # scene.draw_molecule(sfo2.molecule)


    #         scene.draw_isosurface(S, -(.003), opacity=.25, color=[255/255, 105/255, 180/255])
    #         scene.draw_isosurface(S,  (.003), opacity=.25, color=[137/255, 243/255, 54/255])
    #         scene.draw_molecule(sfo1.molecule)
    #         scene.draw_molecule(sfo2.molecule)

    # exit()

    pyfmo.plotting.orbital_diagram.draw_interaction([sfo1, sfo2], 
                                                    [mo1, mo2], 
                                                    [(sfo1, mo1), (sfo1, mo2), (sfo2, mo1), (sfo2, mo2)], 
                                                    orbs, 
                                                    rf'Int. {i+1}: $\Delta E^{{(2)}}_{{ij}} \approx {dE:5.1f}$ kcal/mol (${explanation:5.1%}$% of total)')
    plt.show()

pairs = _get_pauli(orbs, 20)
# mos, sfos, connections  = [], [], []
for i, (dE, explanation, sfo1, sfo2, mo1, mo2) in enumerate(pairs):
    print(dE, explanation, sfo1, sfo2, mo1, mo2)
    # mos.append([mo1, mo2])
    # sfos.append([sfos1, sfos2])
    # connections.append([(sfo1, mo1), (sfo1, mo2), (sfo2, mo1), (sfo2, mo2)])
    pyfmo.plotting.orbital_diagram.draw_interaction([sfo1, sfo2], 
                                                    [mo1, mo2], 
                                                    [(sfo1, mo1), (sfo1, mo2), (sfo2, mo1), (sfo2, mo2)], 
                                                    orbs, 
                                                    rf'Int. {i+1}: $\Delta E^{{(4)}}_{{ij}} \approx {dE:5.1f}$ kcal/mol (${explanation:5.1%}$% of total)')
    plt.show()
exit()


f1, f2 = orbs.fragments
sfos1, sfos2 = orbs.sfos.get_fragment_sfos(f1), orbs.sfos.get_fragment_sfos(f2)
sfos1_occ, sfos2_occ = np.array([sfo.occupied for sfo in sfos1]).reshape(1, -1), np.array([sfo.occupied for sfo in sfos2]).reshape(-1, 1)


E1, E2 = np.array([sfo.energy for sfo in sfos1]).reshape(1, -1), np.array([sfo.energy for sfo in sfos2]).reshape(-1, 1)
dE = abs(E2 - E1)
S = orbs.data.matrices.overlap.total
S = orbs.overlap_matrix(sfos1, sfos2)

occ_virt_mask = np.logical_xor(sfos1_occ, sfos2_occ)
oi = -S**2 / dE * occ_virt_mask
oi = oi / np.sum(oi) * res.properties.energy.orbint.total

rows = []
for i in range(10):
    best = np.unravel_index(argNmax(-oi, i), oi.shape)
    sfo1, sfo2  = sfos1[best[1]], sfos2[best[0]]

    if sfo2.occupied:
        sfo1, sfo2 = sfo2, sfo1

    mos = list(orbs.mos)
    mos = list(sorted(mos, key=lambda mo: mo.energy))
    occmo = np.array([mo.occupied for mo in mos])
    sfo1_contr = np.array([sfo1.mulliken_contribution(mo) for mo in mos])
    sfo2_contr = np.array([sfo2.mulliken_contribution(mo) for mo in mos])

    sfos = list(orbs.sfos)
    sfos = list(sorted(sfos, key=lambda sfo: sfo.energy))

    occ_contrs = abs(sfo1_contr * sfo2_contr) * occmo
    virt_contrs = abs(sfo2_contr * sfo1_contr) * (1-occmo)

    occ_mo = mos[np.argmax(occ_contrs)]
    virt_mo = mos[np.argmax(virt_contrs)]

    pyfmo.plotting.orbital_diagram.draw_interaction([sfo1, sfo2], 
                                                    [occ_mo, virt_mo], 
                                                    [(sfo1, occ_mo), (sfo1, virt_mo), (sfo2, occ_mo), (sfo2, virt_mo)], 
                                                    orbs, 
                                                    f'Explains {oi[best]/res.properties.energy.orbint.total:5.1%} of OI')
    plt.show()

    rows.append([sfos1[best[1]], sfos2[best[0]], occ_mo, virt_mo, f'{oi[best]:5.1f}', f'{oi[best]/res.properties.energy.orbint.total:5.1%}', f'{abs(S[best]):5.1%}', f'{dE[best]:5.2f}'])


tcutility.log.table(rows, ['SFO1', 'SFO2', 'MO1', 'MO2', 'ΔE⁽²⁾', 'Frac.', 'S', 'Δε (eV)'])
plt.figure()
plt.plot(oi.flatten(), c='b')
plt.title('SFO interactions')

log.log()
# print(res.properties.energy.pauli.total)

f1, f2 = orbs.fragments
sfos1, sfos2 = orbs.sfos.get_fragment_sfos(f1), orbs.sfos.get_fragment_sfos(f2)
sfos1_occ, sfos2_occ = np.array([sfo.occupied for sfo in sfos1]).reshape(1, -1), np.array([sfo.occupied for sfo in sfos2]).reshape(-1, 1)

E1, E2 = np.array([sfo.energy for sfo in sfos1]).reshape(1, -1), np.array([sfo.energy for sfo in sfos2]).reshape(-1, 1)
dE = abs(E2 - E1)
S = orbs.data.matrices.overlap.total
S = orbs.overlap_matrix(sfos1, sfos2)

occ_occ_mask = np.logical_and(sfos1_occ, sfos2_occ)
pauli = S**2 * occ_occ_mask
pauli = pauli / np.sum(pauli) * res.properties.energy.pauli.total


rows = []
for i in range(10):
    best = np.unravel_index(argNmax(pauli, i), pauli.shape)
    rows.append([sfos1[best[1]], sfos2[best[0]], f'{pauli[best]:5.1f}', f'{pauli[best]/res.properties.energy.pauli.total:5.1%}', f'{abs(S[best]):5.1%}'])
    # print()

tcutility.log.table(rows, ['SFO1', 'SFO2', 'ΔE⁽⁴⁾', 'Frac.', 'S'])

plt.plot(pauli.flatten(), c='r')
plt.show()


