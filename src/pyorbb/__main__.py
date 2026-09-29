def main():
    import platformdirs
    import os
    # import contextlib
    from datetime import datetime
    import traceback
    # import faulthandler

    # faulthandler.enable()
    
    log_dir = platformdirs.user_log_dir(appname="PyOrbb", appauthor="TheoCheMVU", ensure_exists=True)
    now = str(datetime.now()).replace(" ", "_").replace(":", "-").split(".")[0]
    log_file = os.path.join(log_dir, now + ".txt")
    print("WRITING STDOUT/STDERR TO", log_file)
    # with open(log_file, "w+") as outp:
        # with contextlib.redirect_stdout(outp), contextlib.redirect_stderr(outp):
    try:
        from pyorbb.application import main

        with main.PyOrbbApp():
            ...

    except Exception as e:
        print("".join(traceback.format_exception(type(e), e, e.__traceback__)))


if __name__ == '__main__':
    main()
