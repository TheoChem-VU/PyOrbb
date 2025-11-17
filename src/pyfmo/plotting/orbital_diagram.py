import matplotlib.pyplot as plt
import numpy as np
import pyfmo
import re



def draw_interaction(sfos, mos, connections, 
        title=None, 
        energy_type='energy', 
        connection_colors={}, 
        ax=None, 
        ylim=None, 
        draw_mo_labels=False,
        draw_sfo_labels=True,
        alpha_range=(0.1, 1),
        highlighted_orbitals=None,
        use_darkmode=False,
        font="Aptos",
        fontsize=9,
        **kwargs):

    arrow_length        = kwargs.get('arrow_length', .3 / 4.8280888207)
    arrow_thickness     = kwargs.get('arrow_thickness', .35)
    arrow_width         = kwargs.get('arrow_width', .005)
    arrow_head_width    = kwargs.get('arrow_head_width', .025)
    arrow_head_length   = kwargs.get('arrow_head_length', .1 / 4.8280888207)
    arrow_overhang      = kwargs.get('arrow_overhang', .4)
    arrow_spacing       = kwargs.get('arrow_spacing', .012)

    level_width         = kwargs.get('level_width', .08)
    level_thickness     = kwargs.get('level_thickness', 3)
    highlight_thickness = kwargs.get('highlight_thickness', 2)

    if use_darkmode: 
        highlight_color = '#36B8FF'
        level_color = '#BCBCBC'
        label_color = 'lightgrey'
        axis_label_color = 'white'
        spine_color = 'white'
        
        ax.set_facecolor('#3d3d3d')
        ax.get_figure().patch.set_facecolor('#3d3d3d')
    else:
        highlight_color = '#36B8FF'
        level_color = 'k'
        label_color = 'grey'
        axis_label_color = 'black'
        spine_color = 'k'

    if draw_mo_labels:
        degenerate_mo_threshold = .15
    else:
        degenerate_mo_threshold = .008

    if draw_sfo_labels:
        degenerate_sfo_threshold = .15
    else:
        degenerate_sfo_threshold = .008

    if highlighted_orbitals is None:
        highlighted_orbitals = []

    if ax is None:
        ax = plt.gca()

    if ylim is None:
        try:
            energies = [getattr(orb, energy_type) for orb in list(sfos)] + [orb.energy for orb in list(mos)]
            energy_span = max(energies) - min(energies)
            ax.set_ylim(min(energies) - .1 * energy_span, max(energies) + .1 * energy_span, auto=False)
        except ValueError:
            energy_span = 1
            ax.set_ylim(0, 1, auto=False)
        energy_span *= 1.2
    else:
        energy_span = ylim[1] - ylim[0]
        ax.set_ylim(*ylim, auto=False)


    frags = sorted(set(sfo.fragment for sfo in sfos))
    ax.set_xlim(-1, len(frags), auto=False)
    sep_orbs = {frag: [sfo for sfo in sfos if sfo.fragment == frag] for frag in frags}
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

                E1, E2 = orb.energy, other_orb.energy
                if orb in sfos:
                    E1, E2 = getattr(orb, energy_type), getattr(other_orb, energy_type)

                if isinstance(orb, pyfmo.orbitals.objects.MO):
                    if abs(E1 - E2) < (degenerate_mo_threshold * energy_span):
                        degenerates[-1].append(other_orb)
                else:
                    if abs(E1 - E2) < (degenerate_sfo_threshold * energy_span):
                        degenerates[-1].append(other_orb)

        degenerates = [list(sorted(deg, key=lambda orb: orb.energy)) for deg in degenerates]
        for orb in sep_orbs_:
            orb_degenerate = [deg for deg in degenerates if orb in deg][0]
            deg_idx = orb_degenerate.index(orb) + 1
            deg_degree = len(orb_degenerate) + 1
            poss[orb] = base_pos + 1 / deg_degree * deg_idx

    xtick_pos, xtick_label = [.5], ['Complex']
    for orb in poss:
        if orb not in sfos:
            continue
        if orb.fragment in xtick_label:
            continue

        xtick_pos.append(poss[orb])
        xtick_label.append(orb.fragment)


    ax.set_title(title)
    ax.set_ylabel('Orbital Energy / eV', color=axis_label_color)
    ax.set_xticks(xtick_pos, xtick_label, color=spine_color)
    ax.spines[['top', 'bottom', 'right']].set_visible(False)
    ax.spines['left'].set_color(spine_color)
    ax.yaxis.label.set_color(spine_color)
    ax.tick_params(axis='y', colors=spine_color)

    ax.tick_params('x', labelsize=12, labelcolor=label_color)
    ax.tick_params(bottom = False)
    for orb in poss:
        E = orb.energy
        if orb in sfos:
            E = getattr(orb, energy_type)

        orb_name = pyfmo.generate_label(orb)

        is_MO = orb in mos
        # if our orbital is highlighted we draw an extra plot around it with a different color
        if orb in highlighted_orbitals:
            ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2],
                    [E, E],
                    c=highlight_color,
                    linewidth=level_thickness + highlight_thickness,
                    gid=f'{"MO" if is_MO else "SFO"}_{orb}')

        ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2], 
                [E, E], 
                c=level_color, 
                linewidth=level_thickness, 
                gid=f'{"MO" if is_MO else "SFO"}_{orb}')

        if (is_MO and draw_mo_labels) or (not is_MO and draw_sfo_labels):
            ax.text(poss[orb],
                     E - arrow_length / 1.8 * energy_span,
                     orb_name,
                     ha='center',
                     va='top',
                     size=9,
                     gid=f'{"TEXTMO" if is_MO else "TEXTSFO"}_{orb}',
                     fontname='monospace',
                     color=level_color)

        if not orb.occupied:
            continue

        for spin_part in orb.spin:
            break_on_one = False
            if spin_part == 'A':
                offset_x = -arrow_spacing
                offset_y = -arrow_length / 2 * energy_span
                displacement = arrow_length * energy_span
            elif spin_part == 'B':
                offset_x =  arrow_spacing
                offset_y =  arrow_length / 2 * energy_span
                displacement = -arrow_length * energy_span

            if orb.spin == 'AB' and orb.occupation == 1:
                offset_x = 0
                if orb.spin_pol in (0, 1):
                    offset_y = -arrow_length / 2 * energy_span
                    displacement = arrow_length * energy_span
                elif orb.spin_pol == -1:
                    offset_y = arrow_length / 2 * energy_span
                    displacement = -arrow_length * energy_span
                break_on_one = True

            if orb.spin != 'AB':
                offset_x = 0

            # if orb in highlighted_orbitals:
            #     ax.arrow(poss[orb]+offset_x, 
            #              E+offset_y, 
            #              0, 
            #              displacement, 
            #              width=arrow_width, 
            #              head_width=arrow_head_width, 
            #              head_length=arrow_head_length * energy_span, 
            #              color=highlight_color, 
            #              overhang=arrow_overhang, 
            #              length_includes_head=True,
            #              linewidth=arrow_thickness + highlight_thickness,
            #              gid=f'{"ARROWMO" if is_MO else "ARROWSFO"}_{orb}')

            ax.arrow(poss[orb]+offset_x, 
                     E+offset_y, 
                     0, 
                     displacement, 
                     width=arrow_width, 
                     head_width=arrow_head_width, 
                     head_length=arrow_head_length * energy_span, 
                     color=level_color, 
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

        c = connection_colors.get((sfo, mo), level_color)
        ax.plot([psfo, pmo], [getattr(sfo, energy_type), mo.energy], c=c, linewidth=1.5, alpha=np.clip(sfo.mulliken_contribution(mo), *alpha_range), gid=f'MIX_{sfo} -> {mo}', zorder=-10)
