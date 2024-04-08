import matplotlib.pyplot as plt
import numpy as np
import pyfmo # noqa
from pyfmo.orbitals import matrices # noqa

def plot_property(orbs1, orbs2, prop=None, propargs={}, propkwargs={}, cmap='Greens', title=None, unit=None, use_relname=False, use_indexname=False, scale=None, **kwargs):
    if cmap is None:
        cmap = 'Greens'
        if hasattr(prop, 'cmap'):
            cmap = prop.cmap

    if scale is None:
        scale = 1
        if hasattr(prop, 'scale'):
            scale = prop.scale

    prop_name = ''
    if hasattr(prop, '__name__'):
        prop_name = prop.__name__

    if title is None:
        title = prop_name
        if hasattr(prop, 'title'):
            title = prop.title

    if unit is None:
        if hasattr(prop, 'unit'):
            unit = prop.unit

    if unit is not None:
        unit = '(' + unit + ')'
    else:
        unit = ''
 
    if callable(prop):
        M = prop(orbs1, orbs2, *propargs, **propkwargs)
    else:
        M = prop


    plotname = orbs1[0].spin + ' ' + orbs1[0].kfpath
    plt.figure(figsize=kwargs.get('figsize', (10, 8)), label=f'{prop_name} {plotname}')
    occ_virt_border1 = [i for i in range(1, len(orbs1)) if orbs1[i-1].occupation != orbs1[i].occupation]
    occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
    occ_virt_border2 = [i for i in range(1, len(orbs2)) if orbs2[i-1].occupation != orbs2[i].occupation]
    occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
    plt.imshow(M, origin='lower', cmap=cmap)
    # gridlines
    plt.hlines(y=np.arange(0, len(orbs1))+0.5, xmin=np.full(len(orbs1), -0.5), xmax=np.full(len(orbs1), len(orbs2)-0.5), color="w", linewidth=1.5)
    plt.vlines(x=np.arange(0, len(orbs2))+0.5, ymin=np.full(len(orbs2), -0.5), ymax=np.full(len(orbs2), len(orbs1)-0.5), color="w", linewidth=1.5)
    # occ_virt border lines
    plt.vlines(occ_virt_border2-.5, -.5, len(orbs1)-.5, colors='k', linewidth=2)
    plt.hlines(occ_virt_border1-.5, -.5, len(orbs2)-.5, colors='k', linewidth=2)
    # text inside cells
    for i in range(len(orbs1)):
        for k in range(len(orbs2)):
            val = M[i, k]
            if np.isnan(val):
                continue
            color = 'w' if val > np.nanmax(M) / 2 else 'k'
            plt.gca().text(k, i, f'{val*scale:.2f}', ha="center", va="center", color=color, fontsize=8)

    try:
        psi1 = r'\phi_{' + orbs1[0].fragment_unique_name + r'}'
        psi2 = r'\phi_{' + orbs2[0].fragment_unique_name + r'}'
    except: # noqa
        try:
            psi1 = r'\phi_{' + orbs1[0].moleculename + r'}'
            psi2 = r'\phi_{' + orbs2[0].moleculename + r'}'
        except: # noqa
            psi1 = r'\phi_{' + 'Frag1' + r'}'
            psi2 = r'\phi_{' + 'Frag2' + r'}'

    plt.xlabel('$'+psi2+'$', fontsize=16)
    plt.ylabel('$'+psi1+'$', fontsize=16)
    yticks = range(len(orbs1))
    xticks = range(len(orbs2))
    if use_relname:
        plt.xticks(xticks, [orb.relative_name for orb in orbs2], rotation=90)
        plt.yticks(yticks, [orb.relative_name for orb in orbs1], rotation=0)
    elif use_indexname:
        plt.xticks(xticks, [orb.index_name for orb in orbs2], rotation=90)
        plt.yticks(yticks, [orb.index_name for orb in orbs1], rotation=0)
    else:
        plt.xticks(xticks, [repr(orb) for orb in orbs2], rotation=90)
        plt.yticks(yticks, [repr(orb) for orb in orbs1], rotation=0)
    plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$ ' + unit, fontsize=16)
    plt.tight_layout()

    return



# if __name__ == '__main__':
#     p = '/Users/Tori/PyFMO/test/fixtures/NH3BH3/adf.rkf'
#     orbs = pyfmo.orbitals.Orbitals(p)
#     orbs1 = orbs.sfos[:'Donor(LUMO+4)']
#     orbs2 = orbs.sfos[:'Acceptor(LUMO+4)']
#     heatmap = pyfmo.orbitals.matrices.overlap_mat(orbs1, orbs2)
#     plot_property(orbs1, orbs2, heatmap)
#     plt.show()

