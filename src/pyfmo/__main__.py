import nslog
from pyfmo.application import main


if __name__ == '__main__':
    with main.PyOrbbApp():
        ...


'''
G(32AA) + C(27AA): S = 0.082, dE =  5.94eV, dP1 * dP2 = 0.059 * 0.015 = 0.000885
    S**2/|dE| = 0.001131986532
    OI = 0.000001001808081

G(36AA) + C(10AA): S = 0.237, dE = 24.88eV, dP1 * dP2 = 0.009 * 0.022 = 0.000198
    S**2/|dE| = 0.002257596463
    OI = 0.0000004470041
'''