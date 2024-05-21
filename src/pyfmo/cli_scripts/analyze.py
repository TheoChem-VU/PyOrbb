""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyfmo


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


def main(args: argparse.Namespace):
    orbs = pyfmo.orbitals.Orbitals(args.rkf)
    orbs.write_excel(args.output)

