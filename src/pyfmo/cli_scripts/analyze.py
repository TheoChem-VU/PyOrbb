""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyfmo
import matplotlib.pyplot as plt


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Read orbital information from an ADF calculation and write them to an Excel file and simultaneously create MO diagram."
    subparser = parent_parser.add_parser('analyze', help=desc, description=desc)
    subparser.add_argument("-o", "--output", 
                           type=str, 
                           help="Set the output Excel file to write to.", 
                           default="pyfmo.xlsx")
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to summarize in an Excel file.")
    subparser.add_argument("mos",
                           type=str,
                           help="Set the MOs to view.",
                           default="HOMO, LUMO")


def main(args: argparse.Namespace):
    orbs = pyfmo.orbitals.Orbitals(args.rkf)
    orbs.write_excel(args.output)


    fig = plt.figure(figsize-(12,12))
    mos = orbs.mos[args.mos]
    print(mos)



if __name__ == '__main__':

    rkf = '/Users/torigijzen/PyFMO/test/fixtures/NH3BH3/adf.rkf'

    main(rkf=rkf)



