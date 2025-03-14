""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyfmo
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec



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
        for mix_ in mixer.orbital_interactions(N=10):
            # print(mix_, mix_.nelectrons())
            if mix_.xiaobo_check(oi_s.val):
                mix += mix_

        for mix_ in mixer.pauli_repulsions(N=1):
            if mix_.xiaobo_check(pauli_s.val):
                mix += mix_
        mix.draw_diagram()
        plt.gcf().canvas.draw_idle()

    # plt.subplots()
    # mix = mixer.orbital_interactions(N=1)[0]
    # mix += mixer.pauli_repulsions(N=1)[0]
    # mixes.extend()
    # plt.figure()

    gs = GridSpec(nrows=4, ncols=1, height_ratios=[1, .05, .05, .05])
    oi_ax = plt.gcf().add_subplot(gs[2, 0])
    pauli_ax = plt.gcf().add_subplot(gs[3, 0])

    oi_s = Slider(oi_ax, 'OI', 0.001, mixer.orbital_interactions(N=1)[0].lowest_contribution, valinit=.025)
    pauli_s = Slider(pauli_ax, 'Pauli', 0.001, mixer.pauli_repulsions(N=1)[0].lowest_contribution, valinit=0.2)
    
    oi_s.on_changed(draw_diagram)
    pauli_s.on_changed(draw_diagram)

    plt.gcf().add_subplot(gs[0, 0])
    draw_diagram(mixer)

    # plt.tight_layout()
    plt.show()