import matplotlib.pyplot as plt
import pyorbb
import os
import tcmu
from math import pi

# the systems described are fragment analyses of two HCr fragments
# each rkf file describes the system at different bond distances 
# between the Cr---Cr atoms

# the folder is structured like this:
# ./rkfs
# |-- step1.rkf
# |-- step2.rkf
# |-- ...
# |-- step30.rkf

# Obtain all the adf.rkf files we want to analyse
# sort them by the number in their name
rkfs_folder = './rkfs'
rkf_file_names = [f for f in os.listdir(rkfs_folder) if os.path.isfile(os.path.join(rkfs_folder, f)) and f != '.DS_Store']
rkf_file_names = sorted(rkf_file_names, key=lambda f: int(f.removeprefix('step').removesuffix('.rkf')))

# construct the full paths of the sorted rkf file names
adf_rkf_files = [os.path.join(rkfs_folder, f) for f in rkf_file_names]

# construct an MOTracker object with the obtained adf.rkf files
T = pyorbb.analysis.MOTracker(rkf_files=adf_rkf_files)

# energies of the 12SIGMA and 13SIGMA MOs tracked along the steps
tracked_E_12SIGMA = T.energy('12SIGMA')
tracked_E_13SIGMA = T.energy('13SIGMA')
tracked_E_14SIGMA = T.energy('14SIGMA')
tracked_E_15SIGMA = T.energy('15SIGMA')

# against the bond distance between the two Cr atoms (atoms 1 and 3)
coordinate = T.geometry(1, 3)

plt.figure(figsize=(5,3))
plt.title('Tracked MO energies')
plt.plot(coordinate, tracked_E_12SIGMA, label=r'12$\Sigma$')
plt.plot(coordinate, tracked_E_13SIGMA, label=r'13$\Sigma$')
plt.plot(coordinate, tracked_E_14SIGMA, label=r'14$\Sigma$')
plt.plot(coordinate, tracked_E_15SIGMA, label=r'15$\Sigma$')
plt.xlabel('Cr---Cr distance / Å')
plt.ylabel('MO energy / eV')
plt.legend(frameon=False)
plt.gca().spines[['right', 'top']].set_visible(False)
plt.tight_layout()
plt.savefig('./tracked.png', dpi=500)
plt.close()

# energies of the 12SIGMA and 13SIGMA MOs untracked
untracked_E_12SIGMA = [orbs.mos['12SIGMA'].energy for orbs in T.orbital_objects]
untracked_E_13SIGMA = [orbs.mos['13SIGMA'].energy for orbs in T.orbital_objects]
untracked_E_14SIGMA = [orbs.mos['14SIGMA'].energy for orbs in T.orbital_objects]
untracked_E_15SIGMA = [orbs.mos['15SIGMA'].energy for orbs in T.orbital_objects]
plt.figure(figsize=(5,3))
plt.title('Untracked MO energies')
plt.plot(coordinate, untracked_E_12SIGMA, label=r'12$\Sigma$')
plt.plot(coordinate, untracked_E_13SIGMA, label=r'13$\Sigma$')
plt.plot(coordinate, untracked_E_14SIGMA, label=r'14$\Sigma$')
plt.plot(coordinate, untracked_E_15SIGMA, label=r'15$\Sigma$')
plt.xlabel('Cr---Cr distance / Å')
plt.ylabel('MO energy / eV')
plt.legend(frameon=False)
plt.gca().spines[['right', 'top']].set_visible(False)
plt.tight_layout()
plt.savefig('./untracked.png', dpi=500)
plt.close()


# we generate videos of the tracked MOs to show that they indeed are correctly tracked
transforms = tcmu.geometry.Transform()
transforms.rotate(y=90/180*pi)
pyorbb.analysis.movie.make_orbital_movie('tracked_12SIGMA.mp4', T.track_mo('12SIGMA'), transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('tracked_13SIGMA.mp4', T.track_mo('13SIGMA'), transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('tracked_14SIGMA.mp4', T.track_mo('14SIGMA'), transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('tracked_15SIGMA.mp4', T.track_mo('15SIGMA'), transform=transforms, fps=15)

pyorbb.analysis.movie.make_orbital_movie('untracked_12SIGMA.mp4', [orbs.mos['12SIGMA'] for orbs in T.orbital_objects], transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('untracked_13SIGMA.mp4', [orbs.mos['13SIGMA'] for orbs in T.orbital_objects], transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('untracked_14SIGMA.mp4', [orbs.mos['14SIGMA'] for orbs in T.orbital_objects], transform=transforms, fps=15)
pyorbb.analysis.movie.make_orbital_movie('untracked_15SIGMA.mp4', [orbs.mos['15SIGMA'] for orbs in T.orbital_objects], transform=transforms, fps=15)
