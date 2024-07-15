""" Module containing functions for quickly viewing orbitals via the command line """
import argparse
from pyfmo import orbitals


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and draw the given orbital"
    subparser = parent_parser.add_parser('showmo', help=desc, description=desc)
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to view in an png file.")
    subparser.add_argument("orb",
                           type=str,
                           help="Set the orbital to view")
    subparser.add_argument("-i", "--isovalue",
                           type=float,
                           help="Set the isovalue",
                           default=0.03)
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
    subparser.add_argument("-m", "--material",
                           nargs="?",
                           type=str,
                           help="Set the material",
                           default="shiny")
    subparser.add_argument("-c1", "--color1",
                           nargs="+",
                           type=int,
                           help="Set the color")
    subparser.add_argument("-c2", "--color2",
                           nargs="+",
                           type=int,
                           help="Set the color")

def get_mo_color(mo):
    
    # setting HOMO/LUMO color
    if mo.occupation > 0:
        color1 = (234/255, 51/255, 35/255)
        color2 = (5/255, 23/255, 206/255)
    else:
        color1 = (117/255, 251/255, 253/255)
        color2 = (235/255, 114/255, 46/255)

    return color1, color2


def main(args: argparse.Namespace):
    
    # load orbital and choose an MO to draw
    orbs = orbitals.Orbitals(args.rkf)
    orb = orbs.mos[args.orb]

    color1, color2 = get_mo_color(orb)
    if args.flip:
        color1, color2 = color2, color1

    if args.reverse:
        color1, color2 = color2, color1

    if args.color1 is not None:
        color1 = (args.color1[0]/255, args.color1[1]/255, args.color1[2]/255)

    if args.color2 is not None:
        color2 = (args.color2[0]/255, args.color2[1]/255, args.color2[2]/255)


    orb.draw(color1=color1, color2=color2, isovalue=args.isovalue, gridsize=args.grid, material=args.material) 
    