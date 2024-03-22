import pyfmo
import numpy as np
import matplotlib.pyplot as plt


class Mixing:
    def __init__(self, mos, sfos, z, r=None):
        self.mos = mos
        self.sfos = sfos
        self.z = z
        self.r = r

        self.degree = len(self.mos)
        self.ratio = None if r is None else r/z


    def __str__(self):
        s = 'Mixing('
        s += f'[{", ".join([sfo.make_name(frag_name=True, relative_name=True, spin=True) for sfo in self.sfos])}] -> [{", ".join([mo.relative_name for mo in self.mos])}]'
        s += f', z={self.z:.5f}'
        if self.r:
            s += f', r={self.r:.5f}'
        return s + ')'

    def __contains__(self, other):
        if isinstance(other, pyfmo.orbitals.mo.MO):
            return other in self.mos
        if isinstance(other, pyfmo.orbitals.sfo.SFO):
            return other in self.sfos

    def __gt__(self, other):
        value = self.r or self.z
        other_value = other.r or other.z
        return abs(value) > abs(other_value)


def get_2mixings(orbs, mo):
    frag1 = list(orbs.fragments)[0]
    frag2 = list(orbs.fragments)[1]

    sfos1 = [sfo for sfo in orbs.sfos[frag1] if abs(orbs.mulliken_contribution(sfo, mo)) > 0.01]
    sfos2 = [sfo for sfo in orbs.sfos[frag2] if abs(orbs.mulliken_contribution(sfo, mo)) > 0.01]

    # print(sfos1, sfos2)

    mixings = []
    for i, sfo1 in enumerate(sfos1):
        mul11 = orbs.mulliken_contribution(sfo1, mo)

        for k, sfo2 in enumerate(sfos2):
            mul21 = orbs.mulliken_contribution(sfo2, mo)

            for mo2 in orbs.mos:
                if mo2 is mo:
                    continue

                mul12 = orbs.mulliken_contribution(sfo1, mo2)
                mul22 = orbs.mulliken_contribution(sfo2, mo2)

                z = mul11 * mul12 * mul21 * mul22
                if z > 0.001:
                    mixings.append(Mixing([mo, mo2], [sfo1, sfo2], z))
    return mixings


def get_3mixings(orbs, mo):
    frag1 = list(orbs.fragments)[0]
    frag2 = list(orbs.fragments)[1]

    sfos1 = [sfo for sfo in orbs.sfos[frag1] if abs(orbs.mulliken_contribution(sfo, mo)) > 0.01]
    sfos2 = [sfo for sfo in orbs.sfos[frag2] if abs(orbs.mulliken_contribution(sfo, mo)) > 0.01]

    if len(sfos1) < 2 and len(sfos2) < 2:
        return get_2mixings(orbs, mo)

    mixings2 = get_2mixings(orbs, mo)
    if not mixings2:
        return None

    mixing = max(mixings2)

    mixings = []
    for k in mixing.sfos:
        for r in orbs.mos:
            if r in mixing:
                continue

            Crk = orbs.mulliken_contribution(k, r)

            for l in sfos1 + sfos2:
                if l in mixing:
                    continue

                Crl = orbs.mulliken_contribution(l, r)

                for t in mixing.mos:
                    Ctl = orbs.mulliken_contribution(l, t)
                    Ctk = orbs.mulliken_contribution(k, t)

                    z = Crk * Crl * Ctl * Ctk
                    mixings.append(Mixing(mixing.mos + [r], mixing.sfos + [l], mixing.z, z))

    return mixings


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






