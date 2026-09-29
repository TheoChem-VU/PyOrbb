""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Start the PyOrbb analysis program"
    subparser = parent_parser.add_parser('open', help=desc, description=desc)
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to analyse. If omitted, starts a clean PyOrbb GUI.",
                           nargs='?',
                           default=None)


def main(args: argparse.Namespace):
    import traceback

    try:
        from pyorbb.application import main

        app = main.PyOrbbApp()
        app.__post_init__()

        win = app.add_window()
        if args.rkf is not None:
            win.windows[0].load_analysis(args.rkf)

        app.exec()
        app.shutdown()

    except Exception as e:
        print("".join(traceback.format_exception(type(e), e, e.__traceback__)))

