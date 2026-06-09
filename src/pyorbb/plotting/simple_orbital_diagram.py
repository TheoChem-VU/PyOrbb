import matplotlib.pyplot as plt
import numpy as np
import pyorbb
import itertools


def draw_interaction(sfos, mos, 
        connections=None,
        title=None,
        energy_type='energy',
        connection_colors={},
        ax=None,
        ylim=None,
        draw_mo_labels=False,
        draw_sfo_labels=True,
        alpha_range=(0.05, 0.3),
        merge_non_radical=True):
    arrow_length        = .14 / 4.8280888207
    arrow_thickness     = .05
    arrow_width         = .005
    arrow_head_width    = .02
    arrow_head_length   = .04 / 4.8280888207
    arrow_overhang      = .4
    arrow_spacing       = .012

    mo_level_outline_width = 0.7

    level_width = .055
    level_thickness = 1

    if connections is None:
        connections = [(sfo, mo, abs(sfo.mulliken_contribution(mo))) for sfo, mo in itertools.product(sfos, mos) if abs(sfo.mulliken_contribution(mo)) > 0.1]

    if merge_non_radical:
        # determine which orbitals, if specified, should be merged into one
        merged = [sfo for sfo in sfos if sfo.spin_total_occupation in (0, 1, 2)]
        merged.extend([mo for mo in mos if mo.spin_total_occupation in (0, 1, 2)])

        sfos = [sfo for sfo in sfos if sfo not in merged or sfo.spin in ('A', 'AB')]
        mos = [mo for mo in mos if mo not in merged or mo.spin in ('A', 'AB')]

    merged = list(sfos) + list(mos)
    _connections = []
    for sfo, mo, strength in connections:
        if sfo not in sfos:
            sfo = [_sfo for _sfo in sfos if _sfo.name == sfo.name][0]
        if mo not in mos:
            mo = [_mo for _mo in mos if _mo.name == mo.name][0]

        if any(conn[0] is sfo and conn[1] is mo for conn in _connections):
            continue
        _connections.append((sfo, mo, float(strength)))

    connections = _connections
    # connections = [conn for conn in connections if conn[0] in sfos and conn[1] in mos]

    degenerate_filled_threshold = arrow_length * 1.4
    degenerate_virtual_threshold = arrow_length / 5

    # store all energies here
    energies = {}
    for sfo in sfos:
        energies[sfo] = getattr(sfo, energy_type)
    for mo in mos:
        energies[mo] = mo.energy

    _energies = {}
    for sfo in sfos:
        other_sfos = sfo.spin_match_orbs
        _energies[sfo] = float((getattr(sfo, energy_type) + sum(getattr(other_sfo, energy_type) for other_sfo in other_sfos)) / (len(other_sfos) + 1))

    for mo in mos:
        other_mos = mo.spin_match_orbs
        _energies[mo] = float((mo.energy + sum(other_mo.energy for other_mo in other_mos)) / (len(other_mos) + 1))

    energies.update(_energies)

    if ax is None:
        ax = plt.gca()

    if ylim is None:
        try:
            energy_span = max(energies.values()) - min(energies.values())
            ax.set_ylim(min(energies.values()) - .1 * energy_span, max(energies.values()) + .1 * energy_span, auto=False)
        except ValueError:
            energy_span = 1
            ax.set_ylim(0, 1, auto=False)
    else:
        energy_span = ylim[1] - ylim[0]
        ax.set_ylim(*ylim, auto=False)

    energy_span *= 1.2

    frags = sorted(set(sfo.fragment_unique for sfo in sfos))
    ax.set_xlim(-1, len(frags), auto=False)
    sep_orbs = {frag: [sfo for sfo in sfos if sfo.fragment_unique == frag] for frag in frags}
    sep_orbs['mo'] = mos
    poss = {}

    for typ, sep_orbs_ in sep_orbs.items():
        if typ == 'mo':
            base_pos = 0
        else:
            idx = frags.index(typ)
            base_pos = idx
            if idx == 0:
                base_pos -= 1

        degenerates = []
        for orb in sep_orbs_:
            if any(orb in degenerates_ for degenerates_ in degenerates):
                continue

            degenerates.append([orb])
            for other_orb in sep_orbs_:
                if orb == other_orb:
                    continue

                E1, E2 = energies[orb], energies[other_orb]

                if other_orb.occupied:
                    if abs(E1 - E2) < (degenerate_filled_threshold * energy_span):
                        degenerates[-1].append(other_orb)
                else:
                    if abs(E1 - E2) < (degenerate_virtual_threshold * energy_span):
                        degenerates[-1].append(other_orb)

        degenerates = [list(sorted(deg, key=lambda orb: energies[orb])) for deg in degenerates]
        degenerate_energies = [sum([energies[orb] for orb in deg])/len(deg) for deg in degenerates]
        # remove duplicates in the degenerate lists
        for orb in sep_orbs_:
            deg_members = [(i, deg) for i, deg in enumerate(degenerates) if orb in deg]
            deg_energies = [sum([energies[orb] for orb in deg])/len(deg) for _, deg in deg_members]
            deg_energy_diffs = [abs(energies[orb] - deg_energy) for deg_energy in deg_energies]
            closest_idx = np.argmin(deg_energy_diffs)
            for i, deg_ in enumerate(deg_members):
                if i == closest_idx:
                    continue
                deg_[1].remove(orb)

        degenerates = [list(sorted(deg, key=lambda orb: energies[orb])) for deg in degenerates]
        degenerate_energies = [sum([energies[orb] for orb in deg])/len(deg) for deg in degenerates]

        for orb in sep_orbs_:
            degenerate_idx, orb_degenerate = [(i, deg) for i, deg in enumerate(degenerates) if orb in deg][0]
            deg_idx = orb_degenerate.index(orb) + 1
            deg_degree = len(orb_degenerate) + 1
            poss[orb] = base_pos + 1 / deg_degree * deg_idx
            energies[orb] = degenerate_energies[degenerate_idx]

    xtick_pos, xtick_label = [.5], ['MOs']
    for orb in list(sfos) + list(mos):
        if orb not in sfos:
            continue
        if orb.fragment_unique in xtick_label:
            continue

        xtick_pos.append(poss[orb])
        xtick_label.append(orb.fragment_unique)

    ax.set_title(title)
    ax.set_ylabel('Orbital Energy / eV')
    ax.set_xticks(xtick_pos, xtick_label)
    # ax.set_yticks([], [])
    ax.spines[['top', 'bottom', 'right']].set_visible(False)
    ax.tick_params('x', labelsize=12, labelcolor='grey')
    ax.tick_params(bottom = False)
    for orb in list(sfos) + list(mos):
        E = energies[orb]
        spin = orb.spin

        if merge_non_radical and orb in merged:
            spin = 'AB'


        spin_part = {
            'A': r' $\alpha$',
            'B': r' $\beta$'
        }.get(spin, '')

        if isinstance(orb, pyorbb.orbitals.objects.MO):
            color = 'k'
        else:
            color = '#b3b3b3'

        is_MO = orb in mos
        ax.plot(
            [poss[orb]-level_width/2, poss[orb]+level_width/2], 
            [E, E], 
            c='k', 
            linewidth=level_thickness + mo_level_outline_width, 
            gid=f'{"MO" if is_MO else "SFO"}_{orb}')

        ax.plot(
            [poss[orb]-level_width/2, poss[orb]+level_width/2], 
            [E, E], 
            c=color, 
            linewidth=level_thickness)

        if not orb.occupied:
            continue

        for spin_part in spin:
            break_after_one = False
            if spin == 'AB' and orb.spin_total_occupation == 1:
                offset_x = 0
                offset_y = -arrow_length / 2 * energy_span
                displacement = arrow_length * energy_span
                break_after_one = True
                if orb.spin == 'B':
                    offset_x =  0
                    offset_y =  arrow_length / 2 * energy_span
                    displacement = -arrow_length * energy_span

            elif spin_part == 'A':
                offset_x = -arrow_spacing
                offset_y = -arrow_length / 2 * energy_span
                displacement = arrow_length * energy_span

            elif spin_part == 'B':
                offset_x =  arrow_spacing
                offset_y =  arrow_length / 2 * energy_span
                displacement = -arrow_length * energy_span

            if spin != 'AB':
                offset_x = 0

            ax.arrow(poss[orb]+offset_x, 
                      E+offset_y, 
                      0, 
                      displacement, 
                      width=arrow_width, 
                      head_width=arrow_head_width, 
                      head_length=arrow_head_length * energy_span, 
                      color='k', 
                      overhang=arrow_overhang, 
                      length_includes_head=True,
                      linewidth=arrow_thickness + mo_level_outline_width,
                      gid=f'{"ARROWMO" if is_MO else "ARROWSFO"}_{orb}')

            ax.arrow(poss[orb]+offset_x, 
                      E+offset_y, 
                      0, 
                      displacement, 
                      width=arrow_width, 
                      head_width=arrow_head_width, 
                      head_length=arrow_head_length * energy_span, 
                      color=color, 
                      overhang=arrow_overhang, 
                      length_includes_head=True,
                      linewidth=arrow_thickness,
                      gid=f'{"ARROWMO" if is_MO else "ARROWSFO"}_{orb}')

            if break_after_one:
                break

    for sfo, mo, strength in connections:
        psfo, pmo = poss[sfo], poss[mo]
        if psfo < pmo:
            psfo += level_width/2
            pmo  -= level_width/2
        else:
            psfo -= level_width/2
            pmo  += level_width/2

        c = connection_colors.get((sfo, mo), 'k')
        ax.plot([psfo, pmo], [energies[sfo], energies[mo]], c=c, linewidth=1, alpha=np.clip(strength, *alpha_range), gid=f'MIX_{sfo} -> {mo}', zorder=-10)


if __name__ == '__main__':
    orbs = pyorbb.Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyOrbb/calculations/PyOrb_testing_2022/HomolyticEthane/Ethane.rkf')

    ylim = (-20, 0)
    plt.figure()
    draw_interaction(
        [sfo for sfo in orbs.sfos if ylim[0] < sfo.energy < ylim[1]], 
        [mo for mo in orbs.mos if ylim[0] < mo.energy < ylim[1]])

    plt.figure()
    draw_interaction(
        [sfo for sfo in orbs.sfos if ylim[0] < sfo.energy < ylim[1]], 
        [mo for mo in orbs.mos if ylim[0] < mo.energy < ylim[1]],
        merge_non_radical=False)
    plt.show()
