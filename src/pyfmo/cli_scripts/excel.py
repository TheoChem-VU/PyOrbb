""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyfmo
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.gridspec import GridSpec
from matplotlib import animation



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

    mixer = pyfmo.analysis.mixing.Mixer(orbs, energy_type='energy')

    def draw_diagram(arg):
        plt.cla()
        plt.title('Old way')
        mix = pyfmo.analysis.mixing.Mixing(orbs)
        if oi_b.get_status()[0]:
            for mix_ in mixer.orbital_interactions(N=10):
                # print(mix_, mix_.nelectrons())
                if mix_.xiaobo_check(oi_s.val):
                    mix += mix_

        if pauli_b.get_status()[0]:
            for mix_ in mixer.pauli_repulsions(N=10):
                if mix_.xiaobo_check(pauli_s.val):
                    mix += mix_
        mix.draw_diagram(ax=main_ax)

        props = dict(boxstyle='round', facecolor='white', alpha=1)  # bbox features
        main_ax.txt = main_ax.text(1.03, 0.98, (' '*30 + '\n')*11, transform=main_ax.transAxes, fontsize=8, fontname='monospace', verticalalignment='top', bbox=props)
        plt.gcf().canvas.draw_idle()

    
    # plt.subplots()
    # mix = mixer.orbital_interactions(N=1)[0]
    # mix += mixer.pauli_repulsions(N=1)[0]
    # mixes.extend()
    # plt.figure()
    plt.figure(figsize=[9, 6.5])
    gs = GridSpec(nrows=4, ncols=2, height_ratios=[1, .05, .05, .05], width_ratios=[.1, .7])
    oi_bax = plt.gcf().add_subplot(gs[2, 0])
    pauli_bax = plt.gcf().add_subplot(gs[3, 0])

    oi_bax.axis('off')
    pauli_bax.axis('off')

    oi_sax = plt.gcf().add_subplot(gs[2, 1])
    pauli_sax = plt.gcf().add_subplot(gs[3, 1])

    oi_b = CheckButtons(oi_bax, labels=['Show'], actives=[True])
    pauli_b = CheckButtons(pauli_bax, labels=['Show'], actives=[True])
    
    oi_s = Slider(oi_sax, 'OI', 0.001, mixer.orbital_interactions(N=1)[0].lowest_contribution, valinit=.025, facecolor='g')
    pauli_s = Slider(pauli_sax, 'Pauli', 0.001, mixer.pauli_repulsions(N=1)[0].lowest_contribution, valinit=0.2, facecolor='r')
    
    oi_b.on_clicked(draw_diagram)
    pauli_b.on_clicked(draw_diagram)

    oi_s.on_changed(draw_diagram)
    pauli_s.on_changed(draw_diagram)

    main_ax = plt.gcf().add_subplot(gs[0, :])
    draw_diagram(mixer)

    def on_plot_hover(event):
        # Iterating over each data member plotted
        lines = plt.gca().get_lines()
        lines = sorted(lines, key=lambda line: line.zorder)
        for curve in lines:
            # Searching which data member corresponds to current mouse position
            if curve.contains(event)[0]:
                # main_ax.txt.set_text(' '*40)
                # plt.gca().draw_artist(main_ax.txt)
                # plt.gcf().canvas.blit()
                gid = curve.get_gid()
                s = ''
                if gid.startswith('MO_'):
                    mo = orbs.mos[gid[3:]]
                    s += f'MO'.ljust(31)
                    s += f'\n   {mo}'.ljust(31)
                    s += f'\n   {mo.relative_name}\n'.ljust(31)
                    s += f'\nEnergy     {mo.energy:.2f} eV'.ljust(31)
                    s += f'\nOccupation {mo.occupation}'.ljust(31)
                    s += f'\nSpin       {mo.spin}'.ljust(31)
                    s += f'\nIrrep      {mo.symmetry}'.ljust(31)
                    s += '\n\n\n\n'

                if gid.startswith('SFO_'):
                    sfo = orbs.sfos[gid[4:]]
                    s += f'SFO'.ljust(31)
                    s += f'\n   {sfo}'.ljust(31)
                    s += f'\n   {sfo.relative_name}\n'.ljust(31)
                    s += f'\nFragment   {sfo.fragment_unique}'.ljust(31)
                    s += f'\nEnergy    {sfo.energy: .2f} eV'.ljust(31)
                    s += f'\nPop.      {sfo.gross_population: .3f}'.ljust(31)
                    s += f'\nSpin-pop. {sfo.gross_spin: .3f}'.ljust(31)
                    s += f'\nSpin       {sfo.spin}'.ljust(31)
                    s += f'\nIrrep      {sfo.symmetry}'.ljust(31)
                    s += '\n\n'

                if gid.startswith('MIX_'):
                    sfo = orbs.sfos[gid[4:].split('->')[0].strip()]
                    mo = orbs.mos[gid[4:].split('->')[1].strip()]
                    s += f'SFO'.ljust(31)
                    s += f'\n   {sfo}'.ljust(31)
                    s += f'\n   {sfo.relative_name}\n'.ljust(31)
                    s += f'\nMO'.ljust(31)
                    s += f'\n   {mo}'.ljust(31)
                    s += f'\n   {mo.relative_name}\n'.ljust(31)
                    s += f'\nContr.    {sfo.mulliken_contribution(mo): .2%}'.ljust(30)
                    s += f'\nCoeff.    {sfo.coefficient(mo): .6f}'.ljust(30)
                    s += f'\nSpin       {sfo.spin}'.ljust(30)
                    s += f'\nIrrep      {sfo.symmetry}'.ljust(30)

                main_ax.txt.set_text(s)
                plt.gca().draw_artist(main_ax.txt)
                plt.gcf().canvas.blit()
                # anim.resume()
                # anim.pause()
                # txt.draw(plt.gcf().canvas.renderer)

# plt.tight_layout()

    # txt = plt.text(0, 0, 'test_text')
    # plt.subplots_adjust(left=0.25)
    # anim = animation.ArtistAnimation(plt.gcf(), [(txt,)], interval=1)
    plt.gcf().canvas.mpl_connect('motion_notify_event', on_plot_hover) 

    plt.tight_layout()
    plt.show()
