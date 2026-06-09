import matplotlib.pyplot as plt
import pyorbb
import os


# Obtain all the adf.rkf files we want to analyse
calc_dir = 'HCrCrH'
calc_files = [f for f in os.listdir(calc_dir) if f.startswith('complex')]
calc_files = list(sorted(calc_files, key=lambda f: int(f.split('.')[1].removeprefix('0'))))
adf_rkf_files = [os.path.join(calc_dir, f, 'adf.rkf') for f in calc_files]

# construct an MOTracker object with the obtained adf.rkf files
T = pyorbb.analysis.MOTracker(rkf_files=adf_rkf_files)
plt.figure(figsize=(3,2))

# plot the bond distance between the two chromium atoms (atoms 1 and 3) 
# against the energy of the 12SIGMA and 13SIGMA Mos

plt.plot(T.geometry(1, 3), T.energy('12SIGMA'), label=r'12$\Sigma$')
plt.plot(T.geometry(1, 3), T.energy('13SIGMA'), label=r'13$\Sigma$')
plt.xlabel('Cr---Cr distance / A')
plt.ylabel('MO energy / eV')
plt.gca().spines[['right', 'top']].set_visible(False)
plt.tight_layout()

# plot the bond distance against the 
plt.figure(figsize=(3,2))
plt.plot(T.geometry(1, 3), [orbs.mos['12SIGMA'].energy for orbs in T.orbital_objects], label=r'12$\Sigma$')
plt.plot(T.geometry(1, 3), [orbs.mos['13SIGMA'].energy for orbs in T.orbital_objects], label=r'13$\Sigma$')
plt.xlabel('Cr---Cr distance / A')
plt.ylabel('MO energy / eV')
plt.legend(frameon=False)
plt.gca().spines[['right', 'top']].set_visible(False)
plt.tight_layout()
plt.show()
