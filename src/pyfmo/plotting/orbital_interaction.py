import openpyxl as xl
import pyfmo
from matplotlib import colormaps
from tcutility import ensure_list
import numpy as np


def orbint_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            if sfo1.occupation == sfo2.occupation:
                ret[-1].append(np.NaN)
            else:
                ret[-1].append((sfo1 @ sfo2)**2/abs(sfo1.energy - sfo2.energy))
    return np.array(ret).squeeze()