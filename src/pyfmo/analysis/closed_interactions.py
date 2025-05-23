import pyfmo
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tcutility import ensure_list

font = {'family': 'helvetica',
        'size': 7}

mpl.rc('font', **font)


class Mixing:
    def __init__(self, orbs, mos, sfos, strength):
        self.orbs = orbs
        self.mos = mos
        self.sfos = sfos
        self.strength = strength

    def __str__(self):
        s = f'{self.__class__.__name__}('
        s += f'[{", ".join([sfo.make_name(frag_name=True, relative_name=True, spin=True) for sfo in self.sfos])}]'
        s += ' -> '
        s += f'[{", ".join([mo.relative_name for mo in self.mos])}]'
        s += f', strength={self.strength:.5f})'
        return s

    def __contains__(self, other):
        if isinstance(other, pyfmo.orbitals.mo.MO):
            return other in self.mos
        if isinstance(other, pyfmo.orbitals.sfo.SFO):
            return other in self.sfos

    def __gt__(self, other):
        return abs(self.strength) > abs(other.strength)

    def draw_diagram(self):
        from pyfmo.plotting import orbital_diagram

        frags = list(set(sfo.fragment for sfo in self.sfos))
        frag1, frag2 = frags
        frag1_sfos = [sfo for sfo in self.sfos if sfo.fragment == frags[0]]
        frag2_sfos = [sfo for sfo in self.sfos if sfo.fragment == frags[1]]

        orbital_diagram.diagram(self.mos, self.orbs, frag1_sfos, frag2_sfos)

    def draw_mixing(self, c='red', linewidth=2, label=None):
        if label is None:
            label = str(self)

        frags = list(set(sfo.fragment for sfo in self.sfos))
        frag1, frag2 = frags
        frag1_sfo_idx = [sfo.index for sfo in self.sfos if sfo.fragment == frags[0]]
        frag2_sfo_idx = [sfo.index for sfo in self.sfos if sfo.fragment == frags[1]]
        mo_idx = [mo.index for mo in self.mos]

        sfos1 = self.orbs.sfos[min(frag1_sfo_idx)-2:max(frag1_sfo_idx)+2]
        sfos2 = self.orbs.sfos[min(frag2_sfo_idx)-2:max(frag2_sfo_idx)+2]
        mos = self.orbs.mos[min(mo_idx)-2:max(mo_idx)+2]

        frag1_contr = np.zeros((len(mos), len(sfos1)))
        frag2_contr = np.zeros((len(mos), len(sfos2)))
        for i, mo in enumerate(mos):
            for j, sfo in enumerate(sfos1):
                frag1_contr[i, j] = self.orbs.mulliken_contribution(sfo, mo)
            for j, sfo in enumerate(sfos2):
                frag2_contr[i, j] = self.orbs.mulliken_contribution(sfo, mo)

        fig, axes = plt.subplot_mosaic([[frag1, frag2]])
        plt.suptitle('Closed Interactions')
        im = axes[frag1].imshow(frag1_contr, origin='lower', cmap='Greens', aspect='auto')
        axes[frag1].set_xticks(range(len(sfos1)), [sfo.relative_name for sfo in sfos1], rotation=90)
        axes[frag1].set_yticks(range(len(mos)), [mo.relative_name for mo in mos])
        axes[frag1].set_xlabel(f'{frag1} SFO')
        axes[frag1].set_ylabel('MO')
        offset_grid(axes[frag1])

        im = axes[frag2].imshow(frag2_contr, origin='lower', cmap='Greens', aspect='auto')
        axes[frag2].set_xticks(range(len(sfos2)), [sfo.relative_name for sfo in sfos2], rotation=90)
        axes[frag2].set_yticks([])
        axes[frag2].set_xlabel(f'{frag2} SFO')
        offset_grid(axes[frag2])

        for sfo in self.sfos:
            for mo in self.mos:
                if sfo in sfos1:
                    plt.Rectangle((sfos1.index(sfo)-.5, mos.index(mo)-.5), 1, 1, fill=False, color=c, linewidth=linewidth, zorder=9)
                    axes[frag1].plot((sfos1.index(sfo), len(sfos1)-.5), [mos.index(mo), mos.index(mo)], c=c, label=label, zorder=10)

                else:
                    plt.Rectangle((sfos2.index(sfo)-.5, mos.index(mo)-.5), 1, 1, fill=False, color=c, linewidth=linewidth, zorder=9)
                    axes[frag2].plot((-.5, sfos2.index(sfo)), [mos.index(mo), mos.index(mo)], c=c, label=label, zorder=10)
            
            mo_min = min([mos.index(mo) for mo in self.mos])
            mo_max = max([mos.index(mo) for mo in self.mos])
            if sfo in sfos1:
                axes[frag1].plot((sfos1.index(sfo), sfos1.index(sfo)), [mo_min, mo_max], c=c, linewidth=linewidth, label=label, zorder=10)
            if sfo in sfos2:
                axes[frag2].plot((sfos2.index(sfo), sfos2.index(sfo)), [mo_min, mo_max], c=c, linewidth=linewidth, label=label, zorder=10)

        lims = np.min(np.hstack([frag1_contr, frag2_contr])), np.max(np.hstack([frag1_contr, frag2_contr]))
        im.set_clim(*lims)
        plt.colorbar(im, label='Mulliken Contributions')
        plt.tight_layout()

    @property
    def spin(self):
        spins = list(set(mo.spin for mo in self.mos) & set(sfo.spin for sfo in self.sfos))
        assert len(spins) == 1
        return spins[0]

    @property
    def symmetry(self):
        symmetries = list(set(mo.symmetry for mo in self.mos) & set(sfo.symmetry for sfo in self.sfos))
        assert len(symmetries) == 1
        return symmetries[0]

    @property
    def nocc(self):
        return len([sfo for sfo in self.sfos if sfo.occupied])


class Mixing2(Mixing):
    @property
    def type(self):
        if self.nocc == 0:
            return 'LUMO-LUMO'
        if self.nocc == 1:
            return 'Charge-Transfer'
        if self.nocc == 2:
            return 'Pauli'


class Mixing3(Mixing):
    @property
    def type(self):
        if self.nocc == 0:
            return 'LUMO-LUMO'
        if self.nocc == 1:
            return 'Charge-Transfer'
        if self.nocc == 2:
            return 'Pauli'


def select_mos(orbs, n_mo=10):
    sorted_mos = pyfmo.analysis.orbital_activity.mos(orbs)
    # indices = np.argsort(-activities)
    # sorted_mos = 
    occ_mos = [mo for mo in sorted_mos if mo.occupied][:n_mo]
    virt_mos = [mo for mo in sorted_mos if not mo.occupied][:n_mo]
    return occ_mos + virt_mos


def select_sfos(orbs, n_sfo=10):
    sorted_sfos = pyfmo.analysis.orbital_activity.sfos(orbs)
    # indices = np.argsort(-activities)
    # sorted_sfos = [orbs.sfos.sfos[i] for i in indices]
    ret = {}
    for frag in orbs.fragments:
        frag_sfos = [sfo for sfo in sorted_sfos if sfo.fragment == frag][:n_sfo]
        ret[frag] = frag_sfos
    return ret


def orbint_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            if sfo1.occupation == sfo2.occupation:
                ret[-1].append(np.NaN)
            else:
                ret[-1].append((sfo1 @ sfo2)**2/abs(sfo1.energy - sfo2.energy))
    return np.array(ret).squeeze()

def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


def select_orbs(orbs, n_oi=5, n_pauli=5, n_mo_per_sfo=3):
    # first split the sfos by fragment type
    frag1, frag2 = tuple(orbs.fragments)
    sfos1 = orbs.sfos.get_fragment_sfos(frag1)
    sfos2 = orbs.sfos.get_fragment_sfos(frag2)

    # then select based on OI terms
    oi = orbint_mat(sfos1, sfos2)
    # we set the NaN values to the lowest to prevent some errors
    oi[np.isnan(oi)] = np.nanmin(oi)
    # choose the N highest values for the OIs
    top_N_highest = np.sort(oi.flatten())[-n_oi:]
    # and get their indices
    # here the left column in indices will be for sfos1 
    # and the right column for sfos2
    indices = np.argwhere(np.isin(oi, top_N_highest))
    # get the sfos that correspond to the indices
    oi_sfos = {sfos1[i] for i in indices[:, 0]} | {sfos2[i] for i in indices[:, 1]}

    # for pauli sfos we first get only the occupied sfos
    sfos1_occ = [sfo for sfo in sfos1 if sfo.occupied]
    sfos2_occ = [sfo for sfo in sfos2 if sfo.occupied]

    # we now use the overlap matrix instead of the OI matrix
    S = overlap_mat(sfos1_occ, sfos2_occ)
    # again get the highest N values
    top_N_highest = np.sort(S.flatten())[-n_pauli:]
    # and their indices
    indices = np.argwhere(np.isin(S, top_N_highest))
    # get the sfos that have the biggest overlaps
    pauli_sfos = {sfos1[i] for i in indices[:, 0]} | {sfos2[i] for i in indices[:, 1]}

    # the total SFOs will be the union of the two sets
    selected_sfos = oi_sfos | pauli_sfos

    # with the SFOs selected we can select the MOs
    selected_mos = set()
    # for each SFO we have selected we select the N MOs that it has the highest contributions to
    for sfo in selected_sfos:
        # sort the MOs by contribution of the SFO and select the top N
        mos = set(sorted(orbs.mos.mos, key=lambda mo: -abs(orbs.mulliken_contribution(sfo, mo)))[:n_mo_per_sfo])
        # add the MOs to the set of selected MOs
        selected_mos |= mos

    return list(selected_sfos), list(selected_mos)
    

def get_two_mixing(orbs, n_mo=15, n_sfo=10):
    # split fragment names
    frag1, frag2 = tuple(orbs.fragments)

    # select MOs and SFOs by their activity
    selected_sfos, selected_mos = select_orbs(orbs)
    selected_sfos = {frag: [sfo for sfo in selected_sfos if sfo.fragment == frag] for frag in orbs.fragments}

    mixings = []
    for i, mo1 in enumerate(selected_mos):
        for mo2 in selected_mos[i+1:]:
            for sfo1 in selected_sfos[frag1]:
                for sfo2 in selected_sfos[frag2]:
                    mul11 = orbs.mulliken_contribution(sfo1, mo1)
                    mul21 = orbs.mulliken_contribution(sfo2, mo1)
                    mul12 = orbs.mulliken_contribution(sfo1, mo2)
                    mul22 = orbs.mulliken_contribution(sfo2, mo2)
                    strength = mul11 * mul12 * mul21 * mul22
                    m = Mixing2(orbs, [mo1, mo2], [sfo1, sfo2], strength)
                    mixings.append(m)

    return sorted(mixings, reverse=True)


def get_three_mixing(orbs, N=10, n_sfo=15, n_mo=10):
    # split fragment names
    frag1, frag2 = tuple(orbs.fragments)

    ## select the MOs and SFOs that have the most important contributions
    selected_sfos, selected_mos = select_orbs(orbs)

    mixings = []
    for mix in get_two_mixing(orbs, n_sfo=n_sfo, n_mo=n_mo)[:N]:
        mo1, mo2 = mix.mos
        sfo1, sfo2 = mix.sfos

        for mo3 in selected_mos:
            if mo3 in mix:
                continue

            for sfo3 in selected_sfos:
                if sfo3 in mix:
                    continue

                mul13 = orbs.mulliken_contribution(sfo1, mo3)
                mul23 = orbs.mulliken_contribution(sfo2, mo3)
                mul31 = orbs.mulliken_contribution(sfo3, mo1)
                mul32 = orbs.mulliken_contribution(sfo3, mo2)
                mul33 = orbs.mulliken_contribution(sfo3, mo3)

                strength = mix.strength * mul13 * mul23 * mul31 * mul32 * mul33
                m = Mixing3(orbs, [mo1, mo2, mo3], [sfo1, sfo2, sfo3], strength)
                mixings.append(m)

    return sorted(mixings, reverse=True)


def offset_grid(axis, **kwargs):
    xlim = axis.get_xlim()
    ylim = axis.get_ylim()
    nx = int(xlim[1] - xlim[0])
    ny = int(ylim[1] - ylim[0])

    kwargs.setdefault('color', 'grey')
    kwargs.setdefault('linewidth', .5)

    for x in range(nx):
        axis.plot((x+.5, x+.5), ylim, **kwargs)

    for y in range(ny):
        axis.plot(xlim, (y+.5, y+.5), **kwargs)


if __name__ == '__main__':
    from tcutility import timer

    for _ in range(50):
        # orbs = pyfmo.orbitals2.objects.Orbitals("/Users/yumanhordijk/Downloads/FeCO4CH4.adf.rkf")
        orbs = pyfmo.orbitals2.objects.Orbitals("/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/HydrogenBond/GuanineCytosine.results/adf.rkf")

    # sfo_indices, mo_indices = (orbs.data.matrices.mulliken_contribution.total > 0.1).nonzero()
    # sfos = [orbs.sfos.orbitals[i] for i in sfo_indices]
    # mos = [orbs.mos.orbitals[i] for i in mo_indices]

    # large_contr_pairs = list(zip(sfos, mos))

    # frag = np.array([sfo.fragment for sfo in sfos]).reshape(-1, 1)
    # occ = np.array([sfo.occupation for sfo in sfos]).reshape(-1, 1)
    # occ_match = np.array([sfo.occupied == mo.occupied for sfo, mo in zip(sfos, mos)]).reshape(-1, 1)

    # # for match, pair in zip(occ_match, large_contr_pairs):
    # #     print(pair, match)

    # fragment_mask = (frag != frag.T)
    # occ_mask = (occ != occ.T)
    # occ_match_mask = np.logical_and(occ_match, occ_match.T)
    # mask = np.logical_and(np.logical_and(fragment_mask, occ_mask), occ_match_mask)

    # # C = np.array([sfo.mulliken_contribution(mo) for sfo, mo in large_contr_pairs]).reshape(-1, 1)
    # C = orbs.data.matrices.mulliken_contribution.total[sfo_indices, mo_indices].reshape(-1, 1)
    # S = orbs.data.matrices.overlap.total[sfo_indices, sfo_indices.reshape(-1, 1)]

    # # print(S.shape)
    # # plt.imshow(S)
    # # print(S)
    # # print(sfo_indices)
    # # print(mo_indices)
    # # print(C)
    # # # plt.show()
    # # print(C.shape)

    # # S = np.array()

    # 力 = C @ C.T * S * mask
    # print(力.min(), 力.max())
    # best = np.unravel_index(np.nanargmax(np.abs(力), axis=None), 力.shape)
    # print(large_contr_pairs[best[0]], large_contr_pairs[best[1]])
    # plt.imshow(力)
    # plt.title('POWER')

    # plt.figure()
    # plt.imshow(mask)

    # # [print(pair) for pair in large_contr_pairs]

    # # print(f'{np.sum(orbs.data.matrices.mulliken_contribution.total > 0.001) / orbs.data.matrices.mulliken_contribution.total.size: .1%} taken into account')
    # # plt.figure()
    # # plt.imshow(orbs.data.matrices.mulliken_contribution.total)
    # plt.show()
    # # print(orbs.sfos.orbitals)

    # # # select_orbs(orbs)
    # # mixing = get_two_mixing(orbs)
    # # for mix in mixing:
    # #     print(mix)

    # # mixing = get_three_mixing(orbs)
    # # for mix in mixing:
    # #     print(mix)


    with timer.timer('closed_interactions.calculate_power'):
        sfos = orbs.sfos.orbitals
        frag = np.array([sfo.fragment for sfo in sfos]).reshape(-1, 1)
        occ = np.array([sfo.occupation for sfo in sfos]).reshape(-1, 1)
        eps = np.array([sfo.energy for sfo in sfos]).reshape(-1, 1)
        fragment_mask = (frag != frag.T)
        occ_mask = (occ != occ.T)
        mask = np.logical_and(fragment_mask, occ_mask)

        S = orbs.data.matrices.overlap.total
        oi = S**2 / np.abs(eps - eps.T) * mask
        oi[np.isnan(oi)] = 0

    # plt.imshow(oi)
    # plt.show()
    # print(S[np.abs(S) > 0.001].size / S.size * 100)

    with timer.timer('closed_interactions.select_best'):
        # sfo1_indices, sfo2_indices = (S > 0.001).nonzero()
        best_sfo1, best_sfo2 = np.unravel_index(np.argsort(-np.abs(oi), axis=None), oi.shape)
        print(sfos[best_sfo1[0]], sfos[best_sfo2[0]])

