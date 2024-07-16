""" Module containing functions for quickly viewing orbitals via the command line """
import argparse
from pyfmo import orbitals
import tcviewer


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and draw the given orbital"
    subparser = parent_parser.add_parser('showsfo', help=desc, description=desc)
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
    subparser.add_argument("-o", "--overlap",
                           action="store_true",
                           help="Plot the overlap between f1 and f2")
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
    subparser.add_argument("-sc1", "--sfocolor1",
                           nargs="+",
                           type=int,
                           help="Set the color")
    subparser.add_argument("-sc2", "--sfocolor2",
                           nargs="+",
                           type=int,
                           help="Set the color")
    subparser.add_argument("-oc1", "--overlapcolor1",
                           nargs="+",
                           type=int,
                           help="Set the color")
    subparser.add_argument("-oc2", "--overlapcolor2",
                           nargs="+",
                           type=int,
                           help="Set the color")

def get_sfo_color(sfo):
    # setting HOMO/LUMO color
    if sfo.occupation > 0:
        color1 = (234/255, 51/255, 35/255)
        color2 = (5/255, 23/255, 206/255)
    else:
        color1 = (117/255, 251/255, 253/255)
        color2 = (235/255, 114/255, 46/255)

    return color1, color2


def main(args: argparse.Namespace):
    
    # load orbital and choose an MO to draw
    orbs = orbitals.Orbitals(args.rkf)
    orb = orbs.sfos[args.orb]


    if isinstance(orb, list):
        print(f"There are multiple orbitals associated with '{args.orb}'.")
        print(f"Namely, {orb}.")
        orb = orbs.sfos[input("Enter orbital to visualize: ")]

    color1, color2 = get_sfo_color(orb)
    if args.flip:
        color1, color2 = color2, color1

    if args.reverse:
        color1, color2 = color2, color1

    if args.color1 is not None:
        color1 = (args.color1[0]/255, args.color1[1]/255, args.color1[2]/255)

    if args.color2 is not None:
        color2 = (args.color2[0]/255, args.color2[1]/255, args.color2[2]/255)


    sfo2 = input("Do you want to visualize a second orbital or overlap?: yes/no  ")
    
    if sfo2 == 'no':
        orb.draw(color1=color1, color2=color2, isovalue=args.isovalue, gridsize=args.grid, material=args.material) 
        
    
    if sfo2 == 'yes':

        sfo2 = input("Enter second orbital to visualize: ")
        sfos = orbs.sfos[sfo2]

        if isinstance(sfos, list):
            print(f"There are multiple orbitals associated with '{sfo2}'.")
            print(f"Namely, {sfos}.")
            sfo2 = orbs.sfos[input("Enter orbital to visualize: ")]
        else:
            sfo2 = sfos

        sfo2_color1, sfo2_color2 = get_sfo_color(sfo2)
        if args.flip:
            sfo2_color1, sfo2_color2 = sfo2_color2, sfo2_color1

        if args.sfocolor1 is not None:
            sfo2_color1 = (args.sfocolor1[0]/255, args.sfocolor1[1]/255, args.sfocolor1[2]/255)

        if args.sfocolor2 is not None:
            sfo2_color2 = (args.sfocolor2[0]/255, args.sfocolor2[1]/255, args.sfocolor2[2]/255)


        sfo1 = orb.cube_file(args.grid)
        sfo2 = sfo2.cube_file(args.grid)
        overlap_cub = sfo1.copy()
        overlap_cub.values *= sfo2.values

        with tcviewer.Screen() as scr:

            # visualizes overlap
            if args.overlap:
                overlap_color1 = (139/255, 251/255, 87/255)
                overlap_color2 = (179/255, 47/255, 230/255)

                if args.overlapcolor1 is not None:
                    overlap_color1 = (args.overlapcolor1[0]/255, args.overlapcolor1[1]/255, args.overlapcolor1[2]/255)

                if args.overlapcolor2 is not None:
                    overlap_color2 = (args.overlapcolor2[0]/255, args.overlapcolor2[1]/255, args.overlapcolor2[2]/255)
               
                scr.draw_cub(overlap_cub, color1=overlap_color1, color2 = overlap_color2, isovalue=args.isovalue**2, material=args.material)
            
            else:
                scr.draw_cub(sfo1, color1=color1, color2=color2, isovalue=args.isovalue, material=args.material)
                scr.draw_cub(sfo2, color1=sfo2_color1, color2=sfo2_color2, isovalue=args.isovalue, material=args.material)

