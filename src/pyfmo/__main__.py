import platform
if platform.system() == 'Darwin':
    import nslog
from pyfmo.application import main


if __name__ == '__main__':
    with main.PyOrbbApp():
        ...
