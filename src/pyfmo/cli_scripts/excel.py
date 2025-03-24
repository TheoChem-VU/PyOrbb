""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyfmo
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider, CheckButtons, RadioButtons
from matplotlib.gridspec import GridSpec
from matplotlib.backend_tools import Cursors
from matplotlib import animation
import numpy as np
import os
import tcviewer

font_path = os.path.split(__file__)[0] + '/ibm_plex_mono/IBMPlexMono-Regular.ttf'  # Your font path goes here
mpl.font_manager.fontManager.addfont(font_path)
prop = mpl.font_manager.FontProperties(fname=font_path)

# plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = prop.get_name()


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and write them to an Excel file."
    subparser = parent_parser.add_parser('excel', help=desc, description=desc)
    subparser.add_argument("-o", "--output", 
                           type=str, 
                           help="Set the output Excel file to write to.", 
                           default="pyfmo.xlsx")
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to summarize in an Excel file.")


def main(args: argparse.Namespace):
    orbs = pyfmo.orbitals2.objects.Orbitals(args.rkf)
    orbs.write_excel2(args.output)

    mixers = {etype: pyfmo.analysis.mixing.Mixer(orbs, energy_type=etype) for etype in orbs.sfo_energy_types}

    oi_mixes = {etype: mixer.orbital_interactions(N=100) for etype, mixer in mixers.items()}
    pauli_mixes = {etype: mixer.pauli_repulsions(N=100) for etype, mixer in mixers.items()}

    def update(arg=None):
        draw_diagram(
            ax=main_ax,
            allowed_spins=[s[1:] for s in spin_b.get_checked_labels()],
            allowed_irreps=[i[1:] for i in irrep_b.get_checked_labels()],
            oi_thresh=10**oi_s.val,
            pauli_thresh=pauli_s.val,
            energy_type=orbs.sfo_energy_types[etype_b.index_selected],
            )

    def draw_diagram(ax=None, fig=None, allowed_spins=None, allowed_irreps=None, oi_thresh=None, pauli_thresh=None, ylim=None, energy_type=None):
        if ax is None:
            ax = plt.gca()

        if fig is None:
            fig = plt.gcf()

        ax.clear()
        ax.yaxis.set_major_formatter('{x: 3.0f}')

        global main_mix
        main_mix = pyfmo.analysis.mixing.Mixing(orbs, energy_type=energy_type)
        # print(energy_type)
        if oi_b.get_status()[0]:
            for mix_ in oi_mixes[energy_type]:
                if any(mo.symmetry not in allowed_irreps for mo in mix_.mos):
                    continue
                if any(mo.spin not in allowed_spins for mo in mix_.mos):
                    continue
                if mix_.xiaobo_check(oi_thresh):
                    main_mix += mix_

        if pauli_b.get_status()[0]:
            for mix_ in pauli_mixes[energy_type]:
                if any(mo.symmetry not in allowed_irreps for mo in mix_.mos):
                    continue
                if any(mo.spin not in allowed_spins for mo in mix_.mos):
                    continue
                if mix_.xiaobo_check(pauli_thresh):
                    main_mix += mix_

        main_mix.sanitize()
        main_mix.draw_diagram(ax=main_ax, ylim=ylim)
        sub_mixes = main_mix.split()

        props = dict(edgecolor='white', facecolor='white', alpha=1)  # bbox features
        main_ax.txt = main_ax.text(1.03, 0.98, ' '*37, transform=main_ax.transAxes, fontsize=8, fontname='monospace', verticalalignment='top', bbox=props)
        fig.canvas.draw_idle()


    # plt.subplots()
    # mix = mixer.orbital_interactions(N=1)[0]
    # mix += mixer.pauli_repulsions(N=1)[0]
    # mixes.extend()
    # plt.figure()
    plt.figure(figsize=[9, 6.5])
    gs = GridSpec(nrows=4, ncols=5, height_ratios=[1, .05, .05, .05], width_ratios=[.1, .5, .1, .1, .1])
    oi_bax = plt.gcf().add_subplot(gs[2, 0])
    oi_bax._mouseover_set = set()
    pauli_bax = plt.gcf().add_subplot(gs[3, 0])
    pauli_bax._mouseover_set = set()

    oi_bax.axis('off')
    pauli_bax.axis('off')

    oi_sax = plt.gcf().add_subplot(gs[2, 1])
    pauli_sax = plt.gcf().add_subplot(gs[3, 1])
    oi_sax._mouseover_set = set()
    pauli_sax._mouseover_set = set()

    oi_b = CheckButtons(oi_bax, labels=[' Show'], actives=[True])
    pauli_b = CheckButtons(pauli_bax, labels=[' Show'], actives=[False])
    oi_s_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in oi_mixes.values())
    oi_s = Slider(oi_sax, 'OI', np.log10(0.0001), np.log10(oi_s_max), valinit=np.log10(oi_s_max/1.5), facecolor='g', closedmax=False)
    pauli_s_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in pauli_mixes.values())
    # pauli_s_max = max(mix.xiaobo_value() for mix in pauli_mixes)
    pauli_s = Slider(pauli_sax, 'Pauli', 0.001, pauli_s_max, valinit=pauli_s_max/1.5, facecolor='r')
    
    oi_b.on_clicked(update)
    pauli_b.on_clicked(update)

    oi_s.on_changed(update)
    pauli_s.on_changed(update)

    irrep_ax = plt.gcf().add_subplot(gs[2:4, 2])
    irrep_ax._mouseover_set = set()
    irrep_ax.axis('off')
    irrep_ax.set_title('Irreps')
    irreps = sorted(set([' ' + mo.symmetry for mo in orbs.mos]))
    irrep_b = CheckButtons(irrep_ax, labels=irreps, actives=[True for _ in irreps])
    irrep_b.on_clicked(update)

    spin_ax = plt.gcf().add_subplot(gs[2:4, 3])
    spin_ax._mouseover_set = set()
    spin_ax.axis('off')
    spin_ax.set_title('Spins')
    spins = sorted(set([' ' + mo.spin for mo in orbs.mos]))
    spin_b = CheckButtons(spin_ax, labels=spins, actives=[True for _ in spins])
    spin_b.on_clicked(update)


    etype_ax = plt.gcf().add_subplot(gs[2:4, 4])
    etype_ax._mouseover_set = set()
    etype_ax.axis('off')
    etype_ax.set_title('SFO Energy')
    proper_names = {
        'energy': ' Orbital',
        'site_energy': ' Site',
        'site_energy_SCF0': ' Site (SCF0)',
    }
    etypes = [proper_names[etype]for etype in orbs.sfo_energy_types]
    etype_b = RadioButtons(etype_ax, labels=etypes, active=0)
    etype_b.on_clicked(update)

    main_ax = plt.gcf().add_subplot(gs[0, :])
    update()

    global already_unfaded
    already_unfaded = True

    def _unfade():
        artists = plt.gca().get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)

        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue            

            if not hasattr(artist, 'orig_color'):
                try:
                    artist.orig_color = artist.get_color()
                except:
                    artist.orig_color = artist.get_fc()
            if not hasattr(artist, 'orig_alpha'):
                artist.orig_alpha = artist.get_alpha() or 1

            artist.set_color('white')
            artist.set_alpha(1)
            plt.gca().draw_artist(artist)
            plt.gcf().canvas.blit()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha)
            plt.gca().draw_artist(artist)

            plt.gcf().canvas.blit()

    def _fade_unrelated_ints(orb):
        artists = plt.gca().get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)
        submixes = main_mix.split()
        submix = [submix for submix in submixes if orb in submix.sfos or orb in submix.mos][0]
        faded_artists = []
        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue

            if not hasattr(artist, 'orig_color'):
                try:
                    artist.orig_color = artist.get_color()
                except:
                    artist.orig_color = artist.get_fc()

            if not hasattr(artist, 'orig_alpha'):
                artist.orig_alpha = artist.get_alpha() or 1

            if gid.startswith('MO_'):
                mo = orbs.mos[gid[3:]]
                if mo in submix.mos:
                    continue

            if gid.startswith('SFO_'):
                sfo = orbs.sfos[gid[4:]]
                if sfo in submix.sfos:
                    continue

            if gid.startswith('ARROWMO_'):
                mo = orbs.mos[gid[8:]]
                if mo in submix.mos:
                    continue

            if gid.startswith('ARROWSFO_'):
                sfo = orbs.sfos[gid[9:]]
                if sfo in submix.sfos:
                    continue

            if gid.startswith('TEXTMO_'):
                mo = orbs.mos[gid[7:]]
                if mo in submix.mos:
                    continue

            if gid.startswith('TEXTSFO_'):
                sfo = orbs.sfos[gid[8:]]
                if sfo in submix.sfos:
                    continue

            if gid.startswith('MIX_'):
                sfo = orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = orbs.mos[gid[4:].split('->')[1].strip()]
                if mo in submix.mos:
                    continue
                if sfo in submix.sfos:
                    continue

            faded_artists.append(artist)

        for artist in faded_artists:
            artist.set_color('white')
            artist.set_alpha(1)
            plt.gca().draw_artist(artist)
            plt.gcf().canvas.blit()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha * .1)
            plt.gca().draw_artist(artist)
            plt.gcf().canvas.blit()

        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue
            if artist in faded_artists:
                continue

            artist.set_color('white')
            artist.set_alpha(1)
            plt.gca().draw_artist(artist)
            plt.gcf().canvas.blit()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha)
            plt.gca().draw_artist(artist)

        plt.gcf().canvas.blit()

    def on_plot_hover(event):
        global already_unfaded
        # Iterating over each data member plotted
        lines = plt.gca().get_children()
        lines = sorted(lines, key=lambda line: -line.zorder)
        for curve in lines:
            gid = curve.get_gid()
            if gid is None:
                continue 

            # Searching which data member corresponds to current mouse position
            if not curve.contains(event)[0]:
                continue

            if gid.startswith('MO_') or gid.startswith('SFO_'):
                plt.gcf().canvas.set_cursor(Cursors.HAND)

            s = ''
            if gid.startswith('MO_'):
                mo = orbs.mos[gid[3:]]
                submixes = main_mix.split()
                submix = [submix for submix in submixes if mo in submix.mos][0]
                s += f'MO'.ljust(35)
                s += f'\n   {mo}'
                s += f'\n   {mo.relative_name}\n'
                s += f'\nEnergy     {mo.energy:.2f} eV'
                s += f'\nOccupation {mo.occupation}'
                s += f'\nSpin       {mo.spin}'
                s += f'\nIrrep      {mo.symmetry}'

                s += '\n\nSFO                    Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for sfo in sorted(submix.sfos, key=lambda sfo: -abs(sfo.mulliken_contribution(mo))):
                    s += f'\n{str(sfo):19.19} {sfo.mulliken_contribution(mo): 8.2%} {sfo.coefficient(mo): 7.4f}'

                s += '\n' * (40 - len(s.splitlines()))
                _fade_unrelated_ints(mo)
                already_unfaded = False

                main_ax.txt.set_text(s)
                plt.gca().draw_artist(main_ax.txt)
                plt.gcf().canvas.blit()

                break

            if gid.startswith('SFO_'):
                sfo = orbs.sfos[gid[4:]]
                submixes = main_mix.split()
                submix = [submix for submix in submixes if sfo in submix.sfos][0]
                s += f'SFO'.ljust(35)
                s += f'\n   {sfo}'
                s += f'\n   {sfo.relative_name}\n'
                s += f'\nFragment   {sfo.fragment_unique}'
                s += f'\nEnergy    {getattr(sfo, orbs.sfo_energy_types[etype_b.index_selected]): .2f} eV'
                s += f'\nPop.      {sfo.gross_population: .3f}'
                s += f'\nSpin-pop. {sfo.gross_spin: .3f}'
                s += f'\nSpin       {sfo.spin}'
                s += f'\nIrrep      {sfo.symmetry}'

                s += '\n\nSFO                      S   dE (eV)'
                s += '\n─────────────────── ────── ─────────'
                for sfo2 in sorted(submix.sfos, key=lambda sfo_: -abs(sfo @ sfo_)):
                    if sfo2 == sfo:
                        continue
                    if sfo2.fragment_unique == sfo.fragment_unique:
                        continue
                    s += f'\n{str(sfo2):19.19} {sfo @ sfo2: 5.3f} {abs(sfo.energy - sfo2.energy): 8.2f}'
                
                s += '\n\nMO                     Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for mo in sorted(submix.mos, key=lambda mo: -abs(sfo.mulliken_contribution(mo))):
                    s += f'\n{str(mo):19.19} {sfo.mulliken_contribution(mo): 8.2%} {sfo.coefficient(mo): 7.4f}'

                s += '\n' * (40 - len(s.splitlines()))
                _fade_unrelated_ints(sfo)
                already_unfaded = False

                main_ax.txt.set_text(s)
                plt.gca().draw_artist(main_ax.txt)
                plt.gcf().canvas.blit()

                break

            if gid.startswith('MIX_'):
                sfo = orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = orbs.mos[gid[4:].split('->')[1].strip()]
                connected_sfos = [conn[0] for conn in main_mix.connections if conn[1] == mo and conn[0].fragment_unique != sfo.fragment_unique]
                s += f'SFO'.ljust(35)
                s += f'\n   {sfo}'
                s += f'\n   {sfo.relative_name}\n'
                s += f'\nMO'
                s += f'\n   {mo}'
                s += f'\n   {mo.relative_name}\n'
                s += f'\nContr.    {sfo.mulliken_contribution(mo): .2%}'
                s += f'\nCoeff.    {sfo.coefficient(mo): .6f}'
                s += f'\nSpin       {sfo.spin}'
                s += f'\nIrrep      {sfo.symmetry}'
                s += '\n\nSecond SFO          Bonding?'
                s += '\n─────────────────── ────────'
                for sfo2 in connected_sfos:
                    is_bonding = ((sfo @ sfo2) * sfo.coefficient(mo) * sfo2.coefficient(mo)) >= 0
                    s += f'\n{str(sfo2):19.19} {"   Yes  " if is_bonding else "    No    "}'

                s += '\n' * (40 - len(s.splitlines()))
                _fade_unrelated_ints(mo)
                already_unfaded = False

                main_ax.txt.set_text(s)
                plt.gca().draw_artist(main_ax.txt)
                plt.gcf().canvas.blit()

                break

        else:
            if not already_unfaded:
                already_unfaded = True
                _unfade()

                main_ax.txt.set_text((' '*35 + '\n')*40)
                plt.gca().draw_artist(main_ax.txt)
                plt.gcf().canvas.set_cursor(Cursors.POINTER)
                plt.gcf().canvas.blit()

    global screen
    screen = None
    def on_click(event):
        global screen
        artists = plt.gca().get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)
        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue

            # Searching which data member corresponds to current mouse position
            if not artist.contains(event)[0]:
                continue

            if gid.startswith('MO_'):
                orb = orbs.mos[gid[3:]]

            elif gid.startswith('SFO_'):
                orb = orbs.sfos[gid[4:]]
            else:
                continue

            # orb.draw()
            if screen is None:
                screen = tcviewer.Screen()
                screen.__enter__()
                screen.window.show()

            with screen.add_molscene() as scene:
            # scr.draw_cub(cub, isovalue, material=tcviewer.materials.orbital_shiny)            
                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb.occupied else ([1, .5, 0], [0, 1, 1])
                scene.draw_molecule(orbs.molecule)
                scene.draw_dual_isosurface(orb.cube_file(), colorm=c1, colorp=c2, opacity=.3)
                # scene.draw_isosurface(orb.cube_file(), -0.03, c1, opacity=.3)
                # scene.draw_isosurface(orb.cube_file(),  0.03, c2, opacity=.3)
                scene.draw_text(str(orb))

            screen.exec()

    plt.gcf().canvas.mpl_connect('motion_notify_event', on_plot_hover)
    plt.gcf().canvas.mpl_connect('button_press_event', on_click)

    plt.tight_layout()
    plt.show()
