import matplotlib.pyplot as plt
# from matplotlib import animation
import numpy as np


def draw_interaction(sfos, mos, connections, title=None, energy_type='energy', connection_colors={}, ax=None, ylim=None):
    arrow_length        = .3 / 4.8280888207
    arrow_thickness     = .35
    arrow_width         = .005
    arrow_head_width    = .025
    arrow_head_length   = .1 / 4.8280888207
    arrow_overhang      = .4
    arrow_spacing       = .012

    level_width = .08
    level_thickness = 3
    degenerate_threshold = .08

    if ylim is None:
        try:
            energies = [getattr(orb, energy_type) for orb in list(sfos)] + [orb.energy for orb in list(mos)]
            energy_span = max(energies) - min(energies)
            ax.set_ylim(min(energies) - .1 * energy_span, max(energies) + .1 * energy_span, auto=False)
        except:
            energy_span = 1
            ax.set_ylim(0, 1, auto=False)
    else:
        energy_span = ylim[1] - ylim[0]
        ax.set_ylim(*ylim, auto=False)

    if ax is None:
        ax = plt.gca()

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

                if abs(orb.energy - other_orb.energy) < (degenerate_threshold * energy_span):
                    degenerates[-1].append(other_orb)

        for orb in sep_orbs_:
            orb_degenerate = [deg for deg in degenerates if orb in deg][0]
            deg_idx = orb_degenerate.index(orb) + 1
            deg_degree = len(orb_degenerate) + 1
            poss[orb] = base_pos + 1 / deg_degree * deg_idx

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
    ax.spines[['top', 'bottom']].set_visible(False)
    ax.tick_params('x', labelsize=12, labelcolor='grey')
    ax.tick_params(bottom = False)
    for orb in poss:
        E = orb.energy
        if orb in sfos:
            E = getattr(orb, energy_type)

        spin_part = {
            'A': r' $\alpha$',
            'B': r' $\beta$'
        }.get(orb.spin, '')
        is_MO = orb in mos
        ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2], [E, E], c='k', linewidth=level_thickness, gid=f'{"MO" if is_MO else "SFO"}_{orb}')
        ax.text(poss[orb],
                 E - arrow_length / 1.8 * energy_span,
                 # f'({orb.relative_name.replace("OMO", "").replace("UMO", "")})',
                 f'{orb.name}' + spin_part,
                 ha='center',
                 va='top',
                 size=8,
                 gid=f'{"TEXTMO" if is_MO else "TEXTSFO"}_{orb}',
                 fontname='monospace')

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
        ax.plot([psfo, pmo], [getattr(sfo, energy_type), mo.energy], c=c, linewidth=1, alpha=np.clip(sfo.mulliken_contribution(mo), 0.1, 1), gid=f'MIX_{sfo} -> {mo}', zorder=-10)

    # ax.fill_betweenx(ax.get_ylim(), 0, 1, alpha=.05, facecolor='k')
          