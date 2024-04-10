""" Module containing functions for quickly viewing orbitals via the command line """
import argparse
import pyfmo
from tcviewer import Screen, materials


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and draw the given orbital"
    subparser = parent_parser.add_parser('show', help=desc, description=desc)
    subparser.add_argument("-o", "--output", 
                           type=str, 
                           help="Set the output png file to write to.", 
                           default="pyfmo.png")
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to view in an png file.")
    subparser.add_argument("orb",
                           nargs="?",
                           type=str,
                           help="Set the orbital to view",
                           default="HOMO")
    subparser.add_argument("mat",
                           nargs="?",
                           type=str,
                           help="Set the material to view the orbital",
                           default="tcviewer.materials.orbital_matte")



def main(args: argparse.Namespace):
    
    # load orbital and choose an MO to draw
    orbs = pyfmo.orbitals.Orbitals(args.rkf)
    orb = orbs.mos[args.orb]

    # generate a cub file
    cub = orb.generate_orbital()
    
    with Screen() as scr:
        scr.draw_cub(cub, material=args.mat)
