""" Module containing functions for quickly viewing orbitals via the command line """
import argparse
from pyfmo import orbitals
import tcviewer



def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and draw the given orbital"
    subparser = parent_parser.add_parser('showsfo', help=desc, description=desc)
    subparser.add_argument("-o", "--output", 
                           type=str, 
                           help="Set the output png file to write to.", 
                           default="pyfmo.png")
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
                           type=str,
                           help="Flips the colors")
    subparser.add_argument("-r", "--reverse",
                           type=str,
                           help="Reverse the colors")
    subparser.add_argument("-o", "--overlap",
                           type=str,
                           help="Plot the overlap between f1 and f2")



def main(args: argparse.Namespace):
    
    # load orbital and choose an MO to draw
    orbs = orbitals.Orbitals(args.rkf)
    sfo1 = orbs.sfos[args.sfo1]
    sfo2 = orbs.sfos[args.sfo2]

    if sfo1 != None:
        sfo1[0].draw(isovalue=args.isovalue, gridsize=args.grid) 
    if sfo2 != None:
        sfo2[1].draw(isovalue=args.isovalue, gridsize=args.grid) 

