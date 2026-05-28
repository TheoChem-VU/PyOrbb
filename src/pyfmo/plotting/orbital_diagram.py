import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.bezier import (
    NonIntersectingPathException, get_cos_sin, get_intersection,
    get_parallels, inside_circle, make_wedged_bezier2,
    split_bezier_intersecting_with_closedpath, split_path_inout)
import numpy as np
import pyfmo
import re


def anchored_text(ax, x, y, text, offset_axes=(0.0, 0.1), **kwargs):
    """
    Draw text that is anchored to a point in data space with specified axes offsets.

    Returns the Artist instance.
    """
    fig = ax.figure

    # initial dummy head (will be updated immediately)
    T = ax.text(x, y, text, **kwargs)
    def update_artist(event=None):
        # 1) get tail position in display coords
        anchor_disp = ax.transData.transform((x, y))  # display (pixel) coords

        # 2) convert tail display to axes coords
        # anchor_axes = ax.transAxes.inverted().transform(anchor_disp)  # (0..1, 0..1)
        # update the positions of the arrow
        xyB = anchor_disp[0] + offset_axes[0]/2 * fig.dpi, anchor_disp[1] + offset_axes[1]/2 * fig.dpi
        # we have to set them like this, otherwise it does not work
        # new_disp = ax.transAxes.transform(xyB)
        new_data = ax.transData.inverted().transform(xyB)
        T.set_position(new_data)

    # Connect updates when y-limits change or figure is resized
    ax.callbacks.connect('ylim_changed', update_artist)
    fig.canvas.mpl_connect('resize_event', update_artist)

    # Do an initial update to set correct head position now
    update_artist()

    return T


def arrow_tail_with_axes_offset(ax, anchor, displacement_axes=(0.0, 0.1), anchor_axes=(0.0, 0.0),
                                arrowstyle="<|-", **kwargs):
    """
    Draw an arrow whose tail is fixed to `anchor` in data coords,
    but whose head is offset from the tail by `displacement_axes` measured in
    axes fraction (dx, dy in [0..1] of the axes width/height).

    Returns the ConnectionPatch instance.
    """
    fig = ax.figure

    # initial dummy head (will be updated immediately)
    patch = mpatches.ConnectionPatch(
        anchor, (0, 0), 
        coordsA=ax.transData, coordsB=ax.transData,
        arrowstyle=arrowstyle, shrinkA=0, shrinkB=0, 
        **kwargs
    )
    ax.add_patch(patch)

    def update_patch(event=None):
        # 1) get tail position in display coords
        anchor_disp = ax.transData.transform(anchor)  # display (pixel) coords
        # 2) convert tail display to axes coords
        # anchor_axes = ax.transAxes.inverted().transform(anchor_disp)  # (0..1, 0..1)

        # update the positions of the arrow
        xyA = anchor_disp[0] - displacement_axes[0]/2 * fig.dpi + anchor_axes[0]/2 * fig.dpi, anchor_disp[1] - displacement_axes[1]/2 * fig.dpi + anchor_axes[1]/2 * fig.dpi
        xyB = anchor_disp[0] + displacement_axes[0]/2 * fig.dpi + anchor_axes[0]/2 * fig.dpi, anchor_disp[1] + displacement_axes[1]/2 * fig.dpi + anchor_axes[1]/2 * fig.dpi
        # we have to set them like this, otherwise it does not work
        patch.xy1 = ax.transData.inverted().transform(xyA)
        patch.xy2 = ax.transData.inverted().transform(xyB)

    # Connect updates when y-limits change or figure is resized
    ax.callbacks.connect('ylim_changed', update_patch)
    fig.canvas.mpl_connect('resize_event', update_patch)

    # Do an initial update to set correct head position now
    update_patch()

    return patch


def draw_interaction(sfos, mos, connections, 
        title=None,
        energy_type='energy',
        connection_types={},
        ax=None,
        ylim=None,
        highlighted_orbitals=None,
        use_darkmode=False,
        mo_column_name='Complex',
        xtick_order=None,
        warning_orbs=None,
        **kwargs):
    
    arrow_length        = kwargs.get('arrow_length', 0.25)
    arrow_width         = kwargs.get('arrow_width', .005)
    arrow_head_width    = kwargs.get('arrow_head_width', .025)
    arrow_head_length   = kwargs.get('arrow_head_length', .1 / 4.8280888207)
    arrow_spacing       = kwargs.get('arrow_spacing', .06)
    arrow_color         = kwargs.get('arrow_color', '#000000')

    level_width         = kwargs.get('level_width', .08)
    level_thickness     = kwargs.get('level_thickness', 3)
    highlight_thickness = kwargs.get('highlight_thickness', 2)

    font_name           = kwargs.get('font_name', 'monospace')
    font_size           = kwargs.get('font_size', 9)

    draw_mo_labels      = kwargs.get('draw_mo_labels', False)
    draw_sfo_labels     = kwargs.get('draw_sfo_labels', True)

    orb_label_offset    = kwargs.get('orb_label_offset', -.28)

    alpha_range         = kwargs.get('alpha_range', (0.1, 1))

    degeneracy_threshold = kwargs.get('degeneracy_threshold', 0.15)


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

    frags = []
    for sfo in sfos:
        if sfo.fragment not in frags:
            frags.append(sfo.fragment)

    if xtick_order is None:
        xtick_order = {mo_column_name: .5}
        for i, frag in enumerate(frags):
            if i == 0:
                xtick_order[frag] = -.5
            else:
                xtick_order[frag] = i + .5

    ax.set_xlim(-1, len(frags), auto=False)
    sep_orbs = {frag: [sfo for sfo in sfos if sfo.fragment == frag] for frag in frags}
    sep_orbs[mo_column_name] = mos
    poss = {}

    for typ, sep_orbs_ in sep_orbs.items():
        base_pos = xtick_order[typ] - .5
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

                if abs(E1 - E2) < (degeneracy_threshold * energy_span):
                    degenerates[-1].append(other_orb)

        degenerates = [list(sorted(deg, key=lambda orb: orb.energy)) for deg in degenerates]
        for orb in sep_orbs_:
            orb_degenerate = [deg for deg in degenerates if orb in deg][0]
            deg_idx = orb_degenerate.index(orb) + 1
            deg_degree = len(orb_degenerate) + 1
            poss[orb] = base_pos + 1 / deg_degree * deg_idx

    ax.set_title(title)
    ax.set_ylabel('Orbital Energy / eV', color=axis_label_color)
    ax.set_xticks(list(xtick_order.values()), list(xtick_order.keys()), color=spine_color)
    for i, artist in enumerate(ax.get_xticklabels()):
        if list(xtick_order.keys())[i] == mo_column_name:
            artist.is_MO = True
        else:
            artist.is_MO = False
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

        is_MO = orb in mos

        orb_name = pyfmo.generate_label(orb)
        orb_index = orb.parent.orbitals.index(orb)

        # if our orbital is highlighted we draw an extra plot around it with a different color
        if orb in highlighted_orbitals:
            ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2],
                    [E, E],
                    c=highlight_color,
                    linewidth=level_thickness + highlight_thickness,
                    gid=f'{"MO" if is_MO else "SFO"}_{orb_index}')

        ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2], 
                [E, E], 
                c=level_color, 
                linewidth=level_thickness, 
                gid=f'{"MO" if is_MO else "SFO"}_{orb_index}')

        if (is_MO and draw_mo_labels) or (not is_MO and draw_sfo_labels):
            if warning_orbs is not None:
                for w in warning_orbs:
                    if w[0] is not orb:
                        continue
                    anchored_text(ax,
                             poss[orb],
                             E,
                             '⚠',
                             # [0, -arrow_length / 1.8 * energy_span],
                             [orb_label_offset * 1.8, 0],
                             ha='center',
                             va='center',
                             size=font_size*2.5,
                             clip_on=True,
                             gid=w[1],
                             fontname=font_name,
                             color='r')

            anchored_text(ax,
                     poss[orb],
                     E,
                     orb_name,
                     # [0, -arrow_length / 1.8 * energy_span],
                     [0, orb_label_offset],
                     ha='center',
                     va='top',
                     size=font_size,
                     clip_on=True,
                     gid=f'{"TEXTMO" if is_MO else "TEXTSFO"}_{orb_index}',
                     fontname=font_name,
                     color=level_color)

        if not orb.occupied:
            continue

        for spin_part in orb.spin:
            break_on_one = False
            if spin_part == 'A':
                offset_x = -arrow_spacing
                displacement = arrow_length
            elif spin_part == 'B':
                offset_x =  arrow_spacing
                displacement = -arrow_length

            if orb.spin == 'AB' and orb.occupation == 1:
                offset_x = 0
                if orb.spin_pol in (0, 1):
                    displacement = arrow_length
                elif orb.spin_pol == -1:
                    displacement = -arrow_length
                break_on_one = True

            if orb.spin != 'AB':
                offset_x = 0

            anchor = (poss[orb] + offset_x, E)
            style = mpatches.ArrowStyle.CurveB(
                head_length=arrow_head_length, 
                head_width=arrow_head_width
                )
            ax.add_patch(
                arrow_tail_with_axes_offset(
                    ax,
                    (poss[orb], E),
                    # anchor,
                    [0, displacement],
                    [offset_x, 0],
                    arrowstyle=style,
                    clip_on=True,
                    gid=f'{"ARROWMO" if is_MO else "ARROWSFO"}_{orb_index}',
                    color=arrow_color
                    )
                )

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

        sfo_index = sfo.parent.orbitals.index(sfo)
        mo_index = mo.parent.orbitals.index(mo)
        typ = connection_types.get((sfo, mo), 'Multiple')
        if typ == 'Multiple':
            c = kwargs.get('Multiple_color', '#000000')
        elif typ == 'PR':
            c = kwargs.get('PR_color', '#FF0000')
        elif typ == 'OI':
            c = kwargs.get('OI_color', '#00FF00')
        elif typ == 'Sanitization':
            c = kwargs.get('Sanitization_color', '#b37fb9')

        ax.plot([psfo, pmo], [getattr(sfo, energy_type), mo.energy], c=c, linewidth=1.5, alpha=np.clip(sfo.mulliken_contribution(mo), *alpha_range), gid=f'MIX_{sfo_index} -> {mo_index}', zorder=-10)
