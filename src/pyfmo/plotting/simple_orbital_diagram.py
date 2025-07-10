import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
import numpy as np
import pyfmo


def draw_interaction(sfos, mos, connections,
        title=None,
        energy_type='energy',
        connection_colors={},
        ax=None,
        ylim=None,
        draw_mo_labels=False,
        draw_sfo_labels=True,
        alpha_range=(0.05, 0.3)):
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

    degenerate_filled_threshold = arrow_length * 1.4
    degenerate_virtual_threshold = arrow_length / 5

    # store all energies here
    energies = {}
    for sfo in sfos:
        energies[sfo] = getattr(sfo, energy_type)
    for mo in mos:
        energies[mo] = mo.energy

    if ax is None:
        ax = plt.gca()

    if ylim is None:
        try:
            energy_span = max(energies.values()) - min(energies.values())
            ax.set_ylim(min(energies.values()) - .1 * energy_span, max(energies.values()) + .1 * energy_span, auto=False)
        except:
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
    for orb in poss:
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
    for orb in poss:
        E = energies[orb]
        
        spin_part = {
            'A': r' $\alpha$',
            'B': r' $\beta$'
        }.get(orb.spin, '')

        if isinstance(orb, pyfmo.orbitals.objects.MO):
            orb_name = f'{orb.name}{spin_part}'
            color = '#b3b3b3'
            color = 'k'
        else:
            color = 'k'
            color = '#b3b3b3'
            if orb.spin == 'AB':
                orb_name = orb.name
            else:
                orb_name = f'{orb.name[:-2]}{spin_part}'

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

        for spin_part in orb.spin:
            break_on_one = False
            if orb.spin == 'AB' and orb.occupation == 1:
                offset_x = 0
                offset_y = -arrow_length / 2 * energy_span
                displacement = arrow_length * energy_span
                break_on_one = True

            elif spin_part == 'A':
                offset_x = -arrow_spacing
                offset_y = -arrow_length / 2 * energy_span
                displacement = arrow_length * energy_span
            elif spin_part == 'B':
                offset_x =  arrow_spacing
                offset_y =  arrow_length / 2 * energy_span
                displacement = -arrow_length * energy_span

            if orb.spin != 'AB':
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
            if break_on_one:
                break

    for sfo, mo in connections:
        psfo, pmo = poss[sfo], poss[mo]
        if psfo < pmo:
            psfo += level_width/2
            pmo  -= level_width/2
        else:
            psfo -= level_width/2
            pmo  += level_width/2

        c = connection_colors.get((sfo, mo), 'k')
        ax.plot([psfo, pmo], [energies[sfo], energies[mo]], c=c, linewidth=1, alpha=np.clip(sfo.mulliken_contribution(mo), *alpha_range), gid=f'MIX_{sfo} -> {mo}', zorder=-10)

    # ax.fill_betweenx(ax.get_ylim(), 0, 1, alpha=.05, facecolor='k')
