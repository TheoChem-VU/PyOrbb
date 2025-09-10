import matplotlib.pyplot as plt
import numpy as np
import pyfmo
import re

IRREP_TRANSLATION = {
    "A1": r"A$_1$",
    "A2": r"A$_2$",
    "E1:1": r"A$_1^1$",
    "E1:2": r"A$_1^2$",
    "AA": r"A$^{\prime}$",
    "AAA": r"A$^{\prime\prime}$",
    "A.g": r"A$_g$",
    "A.u": r"A$_u$",
    "B1": r"B$_1$",
    "B2": r"B$_2$",
    "T1": r"T$_1$",
    "T2": r"T$_2$",
    "A1.g": r"A$_\mathrm{1g}$",
    "A1.u": r"A$_\mathrm{1u}$",
    "A2.g": r"A$_\mathrm{2g}$",
    "A2.u": r"A$_\mathrm{2g}$",
    "E.g": r"E$_\mathrm{g}$",
    "E.u":r"E$_\mathrm{u}$",
    "T1.g": r"T$_\mathrm{1g}$",
    "T1.u": r"T$_\mathrm{1u}$",
    "T2.g": r"T$_\mathrm{2g}$", 
    "T2.u": r"T$_\mathrm{2u}$",
    
    "SIGMA": r"$\Sigma$",
    "PI": r"$\Pi$",
    "DELTA": r"$\Delta$",
    "PHI": r"$\Phi$",

    "SIGMA.g": r"$\Sigma_\mathrm{g}$",
    "PI.g": r"$\Pi_g$",
    "DELTA.g": r"$\Delta_\mathrm{g}$",
    "PHI.g": r"$\Phi_g$",
    "SIGMA.u": r"$\Sigma_\mathrm{u}$",
    "PI.u": r"$\Pi_u$",
    "DELTA.u": r"$\Delta_\mathrm{u}$",
    "PHI.u": r"$\Phi_\mathrm{u}$",

    "S": "s",
    "P:x": r"$p_x",
    "P:y": r"$p_y",
    "P:z": r"$p_z",
    "D:xy": r"$d_{xy}$",
    "D:xz": r"$d_{xz}$",
    "D:yz": r"$d_{yz}$",
    "D:z2": r"$d_{z^2}$",
    "D:x2-y2": r"$d_{x^2-y^2}$",
    "F:x": r"$f_{x}$",
    "F:y": r"$f_{y}$",
    "F:z": r"$f_{z}$",
    "F": "$f$",
    "F:xyz": r"$f_{xyz}$",
    "F:z2x": r"$f_{z^2x}$",
    "F:z2y": r"$f_{z^2y}$",
    "F:z3": r"$f_{z^2}$",
    
}

def translate_label(symm_label: str) -> str:
    print(symm_label)
    if symm_label in IRREP_TRANSLATION:
        return IRREP_TRANSLATION[symm_label]
        
    if ":" in symm_label:
        symm, dimension = symm_label.split(":", 1)
        
        if re.search(r"\d+", dimension): 
            dimension = re.sub(r"([xyz])(\d+)", r"\1^{\2}", dimension)
        if symm in IRREP_TRANSLATION:
            if "." in symm:
                symm, parity = symm.split(".", 1)
                return f"{IRREP_TRANSLATION[symm]}$_\mathrm{{{parity}}},\!_{{{dimension}}}$"
            else:
                return f"{IRREP_TRANSLATION[symm]}$_{{{dimension}}}$"

    return symm_label


def draw_interaction(sfos, mos, connections, 
        title=None, 
        energy_type='energy', 
        connection_colors={}, 
        ax=None, 
        ylim=None, 
        draw_mo_labels=False,
        draw_sfo_labels=True,
        alpha_range=(0.1, 1)):
    arrow_length        = .3 / 4.8280888207
    arrow_thickness     = .35
    arrow_width         = .005
    arrow_head_width    = .025
    arrow_head_length   = .1 / 4.8280888207
    arrow_overhang      = .4
    arrow_spacing       = .012

    level_width = .08
    level_thickness = 3
    if draw_mo_labels:
        degenerate_mo_threshold = .08
    else:
        degenerate_mo_threshold = .008

    if draw_sfo_labels:
        degenerate_sfo_threshold = .08
    else:
        degenerate_sfo_threshold = .008

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
    ax.spines[['top', 'bottom', 'right']].set_visible(False)
    ax.tick_params('x', labelsize=12, labelcolor='grey')
    ax.tick_params(bottom = False)
    for orb in poss:
        E = orb.energy
        if orb in sfos:
            # print(orb, E, energy_type)
            E = getattr(orb, energy_type)

        spin_part = {
            'A': r'$\alpha$',
            'B': r'$\beta$'
        }.get(orb.spin, '')

        if isinstance(orb, pyfmo.orbitals.objects.MO):
            orb_name = f'{orb.name}{spin_part}'
            orb_name = orb_name.replace(orb.symmetry, translate_label(orb.symmetry))
        else:
            if orb.spin == 'AB':
                orb_name = orb.name
            else:
                orb_name = f'{orb.name}{spin_part}'
            orb_name = orb_name.replace(orb.subspecies, translate_label(orb.subspecies))

        # if orb.symmetry in IRREP_TRANSLATION:
        #     orb_name = orb_name.replace(orb.symmetry, IRREP_TRANSLATION[orb.symmetry])


        is_MO = orb in mos
        ax.plot([poss[orb]-level_width/2, poss[orb]+level_width/2], [E, E], c='k', linewidth=level_thickness, gid=f'{"MO" if is_MO else "SFO"}_{orb}')

        if (is_MO and draw_mo_labels) or (not is_MO and draw_sfo_labels):
            ax.text(poss[orb],
                     E - arrow_length / 1.8 * energy_span,
                     # f'({orb.relative_name.replace("OMO", "").replace("UMO", "")})',
                     orb_name,
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
                if orb.spin_pol in (0, 1):
                    offset_y = -arrow_length / 2 * energy_span
                    displacement = arrow_length * energy_span
                elif orb.spin_pol == -1:
                    offset_y = arrow_length / 2 * energy_span
                    displacement = -arrow_length * energy_span
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
        ax.plot([psfo, pmo], [getattr(sfo, energy_type), mo.energy], c=c, linewidth=1, alpha=np.clip(sfo.mulliken_contribution(mo), *alpha_range), gid=f'MIX_{sfo} -> {mo}', zorder=-10)
