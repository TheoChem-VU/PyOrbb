""" Module containing functions for quickly viewing orbitals via the command line """
import argparse
from pyfmo import orbitals
from tcutility import ensure_list #noqa


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and draw the given orbital"
    subparser = parent_parser.add_parser('showsfo', help=desc, description=desc)
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to view in an png file.")
    subparser.add_argument("-s1","--sfo1",
                           type=str,
                           help="Set the orbital to view")
    subparser.add_argument("-s2","--sfo2",
                           type=str,
                           help="Set the orbital to view")
    subparser.add_argument("-i", "--isovalue",
                           type=float,
                           help="Set the isovalue",
                           default="0.03")
    subparser.add_argument("-g", "--grid",
                           type=str,
                           help="Set the gridsize",
                           default="medium")
    subparser.add_argument("-f", "--flip",
                           action="store_true",
                           help="Flips the colors")
    subparser.add_argument("-r", "--reverse",
                           action="store_true",
                           help="Reverse the colors")
    subparser.add_argument("-o", "--overlap",
                           action="store_true",
                           help="Plot the overlap between f1 and f2")
    subparser.add_argument("-m", "material",
                           type=str,
                           help="Set the material",
                           default="shiny")

def main(args: argparse.Namespace):
    
    # load orbital and choose an MO to draw
    orbs = orbitals.Orbitals(args.rkf)
    sfo1 = orbs.sfos[args.sfo1]
    sfo2 = orbs.sfos[args.sfo2]

    # plot the given fragment orbitals
    if sfo1 is not None and args.overlap is False:
        sfo1[0].draw(isovalue=args.isovalue, gridsize=args.grid, material=args.material) 
    if sfo2 is not None and args.overlap is False:
        sfo2[1].draw(isovalue=args.isovalue, gridsize=args.grid, material=args.material)

    # plot the overlap
    # if args.overlap == True:

    #     overlap = []

    #     for sfo1[0] in ensure_list(sfo1):
    #         overlap.append([])
    #         for sfo2[1] in ensure_list(sfo2):
    #             overlap[-1].append(abs(sfo1[0] @ sfo2[1]))

        # o1 = sfo1[0]
        # o2 = sfo2[1]

        # # overlap = abs(o1 @ o2)
        # overlap.draw(isovalue=0.009, gridsize=args.grid, material=args.material)









