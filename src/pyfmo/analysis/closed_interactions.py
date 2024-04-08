import pyfmo
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from time import perf_counter

font = {'family': 'helvetica',
        'size': 7}

mpl.rc('font', **font)


class Mixing:
    def __init__(self, mos, sfos, strength):
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

        orbital_diagram.diagram(self.mos, orbs, frag1_sfos, frag2_sfos)

    def draw_mixing(self, c='red', linewidth=2, label=None):
        if label is None:
            label = str(self)

        frags = list(set(sfo.fragment for sfo in self.sfos))
        frag1, frag2 = frags
        frag1_sfo_idx = [sfo.index for sfo in self.sfos if sfo.fragment == frags[0]]
        frag2_sfo_idx = [sfo.index for sfo in self.sfos if sfo.fragment == frags[1]]
        mo_idx = [mo.index for mo in self.mos]

        sfos1 = orbs.sfos[min(frag1_sfo_idx)-2:max(frag1_sfo_idx)+2]
        sfos2 = orbs.sfos[min(frag2_sfo_idx)-2:max(frag2_sfo_idx)+2]
        mos = orbs.mos[min(mo_idx)-2:max(mo_idx)+2]

        frag1_contr = np.zeros((len(mos), len(sfos1)))
        frag2_contr = np.zeros((len(mos), len(sfos2)))
        for i, mo in enumerate(mos):
            for j, sfo in enumerate(sfos1):
                frag1_contr[i, j] = orbs.mulliken_contribution(sfo, mo)
            for j, sfo in enumerate(sfos2):
                frag2_contr[i, j] = orbs.mulliken_contribution(sfo, mo)

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
                    rect = plt.Rectangle((sfos1.index(sfo)-.5, mos.index(mo)-.5), 1, 1, fill=False, color=c, linewidth=linewidth, zorder=9)
                    # axes[frag1].add_patch(rect)
                    handle = axes[frag1].plot((sfos1.index(sfo), len(sfos1)-.5), [mos.index(mo), mos.index(mo)], c=c, label=label, zorder=10)

                else:
                    rect = plt.Rectangle((sfos2.index(sfo)-.5, mos.index(mo)-.5), 1, 1, fill=False, color=c, linewidth=linewidth, zorder=9)
                    # axes[frag2].add_patch(rect)
                    handle = axes[frag2].plot((-.5, sfos2.index(sfo)), [mos.index(mo), mos.index(mo)], c=c, label=label, zorder=10)
            
            mo_min = min([mos.index(mo) for mo in self.mos])
            mo_max = max([mos.index(mo) for mo in self.mos])
            if sfo in sfos1:
                handle = axes[frag1].plot((sfos1.index(sfo), sfos1.index(sfo)), [mo_min, mo_max], c=c, linewidth=linewidth, label=label, zorder=10)
            if sfo in sfos2:
                handle = axes[frag2].plot((sfos2.index(sfo), sfos2.index(sfo)), [mo_min, mo_max], c=c, linewidth=linewidth, label=label, zorder=10)
        # return handle
        lims = np.min(np.hstack([frag1_contr, frag2_contr])), np.max(np.hstack([frag1_contr, frag2_contr]))
        im.set_clim(*lims)
        plt.colorbar(im, label='Mulliken Contributions')
        # fig.legend(handles=handle)
        plt.tight_layout()
        # plt.show()


class Mixing2(Mixing):
    @property
    def type(self):
        nocc = len([sfo for sfo in self.sfos if sfo.occupied])
        if nocc == 0:
            return 'LUMO-LUMO'
        if nocc == 1:
            return 'Charge-Transfer'
        if nocc == 2:
            return 'Pauli'


class Mixing3(Mixing):
    @property
    def type(self):
        nocc = len([sfo for sfo in self.sfos if sfo.occupied])
        if nocc == 0:
            return 'LUMO-LUMO'
        if nocc == 1:
            return 'Charge-Transfer'
        if nocc == 2:
            return 'Pauli'


def select_mos(orbs, n_mo=10):
    activities = pyfmo.analysis.selection.mo_activity(orbs)
    indices = np.argsort(-activities)
    occ_mos = [orbs.mos.mos[i] for i in indices if orbs.mos.mos[i].occupied][:n_mo]
    virt_mos = [orbs.mos.mos[i] for i in indices if not orbs.mos.mos[i].occupied][:n_mo]
    return occ_mos + virt_mos


def select_sfos(orbs, n_sfo=10):
    activities = pyfmo.analysis.selection.sfo_activity(orbs)
    indices = np.argsort(-activities)
    sorted_sfos = [orbs.sfos.sfos[i] for i in indices]
    ret = {}
    for frag in orbs.fragments:
        frag_sfos = [sfo for sfo in sorted_sfos if sfo.fragment == frag][:n_sfo]
        ret[frag] = frag_sfos
    return ret


def get_two_mixing(orbs, n_mo=15, n_sfo=10):
    # split fragment names
    frag1, frag2 = tuple(orbs.fragments)

    # select MOs and SFOs by their activity
    selected_sfos = select_sfos(orbs, n_sfo=n_sfo)
    selected_mos = select_mos(orbs, n_mo=n_mo)

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
                    m = Mixing2([mo1, mo2], [sfo1, sfo2], strength)
                    mixings.append(m)

    return sorted(mixings, reverse=True)


def get_three_mixing(orbs, N=10, n_sfo=15, n_mo=10):
    ## select the MOs and SFOs that have the most important contributions
    selected_sfos = select_sfos(orbs, n_sfo=n_sfo)
    selected_sfos = sum(selected_sfos.values(), start=[])
    selected_mos = select_mos(orbs, n_mo=n_mo)

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
                m = Mixing3([mo1, mo2, mo3], [sfo1, sfo2, sfo3], strength)
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
    orbs = pyfmo.orbitals.Orbitals('/Users/yumanhordijk/PhD/TheoCheM_stack/calculations/PyOrb_testing_2022/TransitionState/DielsAlder.results/adf.rkf')
    # orbs = pyfmo.orbitals.Orbitals('/Users/yumanhordijk/PhD/TheoCheM_stack/PyFMO/test/fixtures/NH3BH3/adf.rkf')
    orbs.write_excel('/Users/yumanhordijk/PhD/TheoCheM_stack/calculations/PyOrb_testing_2022/TransitionState/DielsAlder.results/pyorb.xlsx')
    mixing = get_two_mixing(orbs)
    print(len(mixing))
    for mix in sorted(mixing, reverse=True):
        print(mix)

    mixing = get_three_mixing(orbs)
    print(len(mixing))
    for mix in sorted(mixing, reverse=True):
        print(mix)

    exit()

    plt.figure()
    max(get_two_mixing(mos, sfos)).draw_mixing()

    # plt.figure()
    # max(get_two_mixing(mos, sfos)).draw_diagram()

    plt.figure()
    max(three_mixing(mos, sfos)).draw_mixing()

    # plt.figure()
    # max(three_mixing(mos, sfos)).draw_diagram()

    plt.tight_layout()
    plt.show()

    exit()


    highest = np.zeros((25, 25))
    timings = np.zeros((25, 25))
    # timings = []

    for i in range(25):
        for j in range(25):
            start = perf_counter()
            mos = select_mos(orbs, j+3)
            sfos = select_sfos(orbs, i+3)
            mix2 = three_mixing(mos, sfos)
            try:
                # highest.append(max(mix2).z)
                highest[i,j] = max(mix2).z
                print(max(mix2))
            except ValueError:
                highest[i,j] = 0
            print(f'Time taken [{i=}, {j=}]: {perf_counter() - start: .2f}s')
            timings[i,j] = perf_counter() - start

    plt.imshow(highest, origin='lower')
    plt.figure()

    plt.imshow(timings, origin='lower')
    plt.show()
