import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
from PIL.ImageDraw import ImageDraw
from PIL import Image
import os

margin = 30

icon_index = 0


def dist_to_nonzero(arr):
    def get_closest_x(i, j):
        # get left
        row = arr[i, :]
        left = 0
        for k, x in enumerate(row[0:j][::-1]):
            if x != 0:
                left = k
                break

        right = 0
        for k, x in enumerate(row[j:]):
            if x != 0:
                right = k
                break

        # print(left, right)
        return min(left, right)

    def get_closest_y(i, j):
        # get left
        row = arr.T[j, :]
        left = 0
        for k, x in enumerate(row[0:i][::-1]):
            if x != 0:
                left = k
                break

        right = 0
        for k, x in enumerate(row[i:]):
            if x != 0:
                right = k
                break

        # print(left, right)
        return min(left, right)


    ret = np.zeros_like(arr)
    for i, row in enumerate(arr):
        for j, v in enumerate(row):
            if v != 0:
                continue
            ret[i, j] = min(get_closest_x(i, j), get_closest_y(i, j))
    return ret




def make_icon( 
        Nsteps=3,
        orbital_exponent=-8,
        border=40,
        border_style=None,
        border_style_settings={},
        bground_style=None,
        cmap='twilight',
        cmap_background=None,
        psi_mode='d'
        ):
    
    resolution = 800
    x, y = np.arange(0, resolution), np.arange(0, resolution)
    X, Y = np.meshgrid(x, y)

    if psi_mode == 'p':
        psi = np.sqrt((X/resolution*2-1)**2 + (Y/resolution*2-1)**2) * ((X/resolution*2-1) + (Y/resolution*2-1)) * np.exp(orbital_exponent*((X/resolution*2-1)**2 + (Y/resolution*2-1)**2))

    if psi_mode == 'd':
        psi = ((X/resolution*2-1) * (Y/resolution*2-1)) * np.exp(orbital_exponent*((X/resolution*2-1)**2 + (Y/resolution*2-1)**2))

    # psi = np.power(psi, 1.5)
    psi = psi / psi.max() * Nsteps
    psi = np.round(psi)

    im = Image.new('RGB', (resolution, resolution), color=0)
    drawer = ImageDraw(im)
    squircle = drawer.rounded_rectangle(xy=[(margin, margin), (resolution-margin, resolution-margin)], radius=int(120/512 * resolution), fill=(255, 255, 255), outline=(1, 1, 1), width=border, corners=None)
    squircle = np.array(im.getdata()).reshape(im.size[0], im.size[1], 3)[:, :, 0]

    psi[squircle == 0] = None
    psi_bg = None

    if bground_style == 'hide':
        psi[psi == 0] = None

    elif bground_style == 'lowest':
        psi[psi == 0] = 0

    elif bground_style == 'highest':
        psi[psi == 0] = Nsteps

    elif bground_style == 'shaded':

        squircle_bg = drawer.rounded_rectangle(xy=[((margin), (margin)), (resolution-(margin), resolution-(margin))], radius=int(120/512 * resolution), fill=(255, 255, 255), outline=(1, 1, 1), width=5, corners=None)
        squircle_bg = np.array(im.getdata()).reshape(im.size[0], im.size[1], 3)[:, :, 0]
        psi_bg = psi.copy()
        # psi_bg[psi != 0] = None
        # psi_bg[psi == 0] = 0
        # psi[psi==0] = None
        # psi[psi==0] = None

        psi_bg[psi!=0] = None
        a = 12
        b = 15
        psi_bg[psi==0] = (np.exp(a*Y[psi==0]/resolution) - 1) / (np.exp(a) - 1) * 1 + (np.exp(b*X[psi==0]/resolution) - 1) / (np.exp(b) - 1) * 1
        # plt.imshow(psi, cmap='Greys')
        # plt.show()

        d = dist_to_nonzero(psi)
        d = 1/d
        d[np.isinf(d)] = 0
        psi_bg += d *0.5

        psi_bg = psi_bg - np.nanmin(psi_bg)

        # plt.imshow(d)
        # plt.show()
        psi[psi == 0] = None

    if border_style == 'wavy':
        psi[squircle == 1] = (border_style_settings.get('f', Nsteps) * (np.sin(X/border_style_settings.get('sinf', 300)) * np.cos(Y/border_style_settings.get('sinf', 140))))[squircle == 1]

    if border_style == 'highest':
        psi[squircle == 1] = Nsteps

    if border_style == 'lowest':
        psi[squircle == 1] = 0

    px = 1/plt.rcParams['figure.dpi']  # pixel in inches
    plt.figure(figsize=(resolution*px, resolution*px))
    plt.imshow(psi, cmap=cmap, vmin=-Nsteps, vmax=Nsteps)
    # plt.show()
    if psi_bg is not None:
        # psi_bg = psi
        # psi_bg[psi != 0] = None
        # psi_bg[psi == 0] = 0
        plt.imshow(psi_bg, cmap='Greys', vmin=0, vmax=1)

    # plt.imshow(psi, cmap='berlin', vmin=-Nsteps, vmax=Nsteps)
    # plt.imshow(psi, cmap='managua', vmin=-Nsteps, vmax=Nsteps)
    # plt.imshow(psi, cmap='bwr', vmin=-Nsteps, vmax=Nsteps)
    plt.axis('off')
    global icon_index
    icon_index += 1
    os.makedirs('icons', exist_ok=True)
    filename = os.path.join('icons', f'icon_{icon_index}.png')
    print('Saved', filename)
    plt.savefig(filename, bbox_inches='tight', transparent=True)
    plt.close()


for psi_mode in ['d', 'p']:
    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
    make_icon(psi_mode=psi_mode, cmap='berlin', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
    make_icon(psi_mode=psi_mode, cmap='bwr', orbital_exponent=-8, border_style='wavy', border=40, bground_style='lowest')
    make_icon(psi_mode=psi_mode, cmap='managua', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
    make_icon(psi_mode=psi_mode, cmap='seismic', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')

    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='wavy', border=30, bground_style='lowest')
    make_icon(psi_mode=psi_mode, cmap='berlin', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
    make_icon(psi_mode=psi_mode, cmap='bwr', orbital_exponent=-8, border_style=None, border=20, bground_style='lowest')

    make_icon(psi_mode=psi_mode, cmap='bwr', bground_style='shaded', orbital_exponent=-7, Nsteps=3, border_style=None, border=20)
    make_icon(psi_mode=psi_mode, cmap='berlin', bground_style='shaded', orbital_exponent=-7, Nsteps=3, border_style=None, border=20)
    make_icon(psi_mode=psi_mode, cmap='seismic', bground_style='shaded', orbital_exponent=-7, Nsteps=3, border_style=None, border=20)
    make_icon(psi_mode=psi_mode, cmap='twilight', bground_style='shaded', orbital_exponent=-7, Nsteps=3, border_style=None, border=20)
    make_icon(psi_mode=psi_mode, cmap='managua', bground_style='shaded', orbital_exponent=-7, Nsteps=3, border_style=None, border=20)

    make_icon(psi_mode=psi_mode, cmap='managua', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
    make_icon(psi_mode=psi_mode, cmap='seismic', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')

    make_icon(psi_mode=psi_mode, cmap='berlin', orbital_exponent=-8, border_style=None, border=20, bground_style='lowest')

    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='highest', border=30, bground_style='lowest')
    make_icon(psi_mode=psi_mode, Nsteps=3, orbital_exponent=-8, border_style='lowest', border=30, bground_style='lowest')
    make_icon(psi_mode=psi_mode, Nsteps=2, orbital_exponent=-7, border_style='lowest', border=30, bground_style='lowest')
    make_icon(psi_mode=psi_mode, Nsteps=20, orbital_exponent=-13, border_style='wavy', border=30, bground_style='lowest')
    make_icon(psi_mode=psi_mode, Nsteps=20, orbital_exponent=-10, border_style=None, border=30, bground_style='lowest')

    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='highest', border=30, bground_style='highest')
    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='lowest', border=30, bground_style='highest')
    make_icon(psi_mode=psi_mode, orbital_exponent=-8, border_style='wavy', border=30, bground_style='highest')
