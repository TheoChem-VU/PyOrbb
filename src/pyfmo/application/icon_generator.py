import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
from PIL.ImageDraw import ImageDraw
from PIL import Image

margin = 30

def make_icon(
		filename, 
		Nsteps=3,
		orbital_exponent=-8,
		border=40,
		border_style=None,
		border_style_settings={},
		bground_style=None,
		cmap='twilight',
		):

	x, y = np.arange(0, 512), np.arange(0, 512)
	X, Y = np.meshgrid(x, y)

	psi = (X/512*2-1) * (Y/512*2-1) * np.exp(orbital_exponent*((X/512*2-1)**2 + (Y/512*2-1)**2))

	# psi = np.power(psi, 1.5)
	psi = psi / psi.max() * Nsteps
	psi = np.round(psi)

	im = Image.new('RGB', (512, 512), color=0)
	drawer = ImageDraw(im)
	squircle = drawer.rounded_rectangle(xy=[(margin, margin), (512-margin, 512-margin)], radius=120, fill=(255, 255, 255), outline=(1, 1, 1), width=border, corners=None)
	squircle = np.array(im.getdata()).reshape(im.size[0], im.size[1], 3)[:, :, 0]

	psi[squircle == 0] = None


	if bground_style == 'hide':
		psi[psi == 0] = None

	if bground_style == 'lowest':
		psi[psi == 0] = 0

	if border_style == 'wavy':
		psi[squircle == 1] = (border_style_settings.get('f', 5) * (np.sin(X/border_style_settings.get('sinf', 200)) * np.cos(Y/border_style_settings.get('sinf', 140))))[squircle == 1]

	plt.imshow(psi, cmap=cmap, vmin=-Nsteps, vmax=Nsteps)
	# plt.imshow(psi, cmap='berlin', vmin=-Nsteps, vmax=Nsteps)
	# plt.imshow(psi, cmap='managua', vmin=-Nsteps, vmax=Nsteps)
	# plt.imshow(psi, cmap='bwr', vmin=-Nsteps, vmax=Nsteps)
	plt.axis('off')
	plt.savefig(filename, bbox_inches='tight', transparent=True)


make_icon('icon.png', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
make_icon('icon2.png', cmap='berlin', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
make_icon('icon3.png', cmap='bwr', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
make_icon('icon4.png', cmap='managua', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')
make_icon('icon9.png', cmap='seismic', orbital_exponent=-8, border_style='wavy', border=40, bground_style='hide')

make_icon('icon5.png', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
make_icon('icon6.png', cmap='berlin', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
make_icon('icon7.png', cmap='bwr', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
make_icon('icon8.png', cmap='managua', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')
make_icon('icon10.png', cmap='seismic', orbital_exponent=-8, border_style='wavy', border=20, bground_style='lowest')

make_icon('icon11.png', cmap='berlin', orbital_exponent=-8, border_style=None, border=20, bground_style='lowest')
