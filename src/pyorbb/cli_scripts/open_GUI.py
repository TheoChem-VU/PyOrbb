""" Module containing functions for quickly submitting geometry optimization jobs via the command line """
import argparse
import pyorbb


def create_subparser(parent_parser: argparse.ArgumentParser):
    desc = "Start the PyOrbb analysis program"
    subparser = parent_parser.add_parser('open', help=desc, description=desc)
    subparser.add_argument("rkf",
                           type=str,
                           help="The path to the `adf.rkf` file to analyse. If omitted, starts a clean PyOrbb GUI.",
                           nargs='?',
                           default=None)


def main(args: argparse.Namespace):
    import platformdirs
    import os
    from datetime import datetime
    import traceback

    log_dir = platformdirs.user_log_dir(appname="PyOrbb", appauthor="TheoCheMVU", ensure_exists=True)
    now = str(datetime.now()).replace(" ", "_").replace(":", "-").split(".")[0]
    log_file = os.path.join(log_dir, now + ".txt")
    with open(log_file, "w+") as outp:
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

