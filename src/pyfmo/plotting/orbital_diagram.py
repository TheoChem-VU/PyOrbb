import matplotlib.pyplot as plt
# from matplotlib import animation
import numpy as np


def draw_interaction(sfos, mos, connections, title=None, energy_type='energy', connection_colors={}, ax=None):
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

    energies = [getattr(orb, energy_type) for orb in list(sfos)] + [orb.energy for orb in list(mos)]
    energy_span = max(energies) - min(energies)

    if ax is None:
        ax = plt.gca()

    ax.set_ylim(min(energies) - .1 * energy_span, max(energies) + .1 * energy_span, auto=False)
    energy_span *= 1.2

    frags = sorted(set(sfo.fragment_unique for sfo in sfos))
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

    # for mo in mos:
    #     mo_degenerate = [deg for deg in degenerates if mo in deg][0]
    #     deg_idx = mo_degenerate.index(mo) + 1
    #     deg_degree = len(mo_degenerate) + 1

    #     poss[mo] = 1 / deg_degree * deg_idx

    # # poss = {mo: .5 for mo in mos}
    # frags = sorted(set(sfo.fragment_unique for sfo in sfos))
    # for sfo in sfos:
    #     idx = frags.index(sfo.fragment_unique)
    #     if idx == 0:
    #         poss[sfo] = idx - .5
    #     else:
    #         poss[sfo] = idx + .5

    # sfos_by_frag = {frag: [sfo for sfo in sfos if sfo.fragment_unique == frag] for frag in frags}
    # degenerates = []
    # for frag, frag_sfos in sfos_by_frag.items():
    #     for sfo in frag_sfos:
    #         if any(sfo in degenerates_ for degenerates_ in degenerates):
    #             continue

    #         degenerates.append([sfo])
    #         for other_sfo in sfos:
    #             if sfo == other_sfo:
    #                 continue

    #             if abs(sfo.energy - other_sfo.energy) < (degenerate_threshold * energy_span):
    #                 degenerates[-1].append(other_sfo)

    # for frag, frag_sfos in sfos_by_frag.items():
    #     for sfo in frag_sfos:
    #         sfo_degenerate = [deg for deg in degenerates if sfo in deg][0]
    #         deg_idx = sfo_degenerate.index(sfo) + 1
    #         deg_degree = len(sfo_degenerate) + 1

    #         poss[sfo] += 1 / deg_degree * deg_idx - .5

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
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
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
        ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2], [E, E], c='k', linewidth=level_thickness, gid=f'{"MO_" if is_MO else "SFO_"}{orb}')
        ax.text(poss[orb],
                 E - arrow_length / 1.8 * energy_span,
                 # f'({orb.relative_name.replace("OMO", "").replace("UMO", "")})',
                 f'{orb.name}' + spin_part,
                 ha='center',
                 va='top',
                 size=6)

        if not orb.occupied:
            continue

        for spin_part in orb.spin:
            if spin_part == 'A':
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
                      linewidth=arrow_thickness)

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

    ax.fill_betweenx((min(energies) - .1 * energy_span, max(energies) + .1 * energy_span), 0, 1, alpha=.05, facecolor='k')
          
    # plt.show()
