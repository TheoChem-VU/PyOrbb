'''
This example shows how to generate a video of the evolution of the 
orbital interaction mechanism of the oxidative addition of a 
palladium atom into the methane C-H bond using PyOrbb.

This example uses ADF fragment calculations along the IRC path of 
the oxidative addition reaction. All data was computed at ZORA-BP86/TZ2P 
with Good numerical quality.
'''

import pyfmo
import os
import matplotlib.pyplot as plt
import numpy as np
import moviepy.editor as mvp


# obtain all adf.rkf files
files = [os.path.join('OxAdd_rkfs', file) for file in os.listdir('OxAdd_rkfs') if file.endswith('.adf.rkf')]

# then load the Orbitals objects
orbs = [pyfmo.Orbitals(file) for file in files]

# we want to also have the key distance, the C-H distance
# we will sort on this distance later on
molecules = [orb.molecule for orb in orbs]
# C and H are atoms 1 and 3 in the molecule
distances = [mol[1].distance_to(mol[3]) for mol in molecules]
# we then get the order in which to generate the frames
distance_order = np.argsort(distances)

# sort the orbs and distances by the order
orbs = [orbs[i] for i in distance_order]
distances = [distances[i] for i in distance_order]

# make the directory to store the frames in
os.makedirs('frames', exist_ok=True)

# for each Orbitals object generate the diagram and save it to a file
# we want to also show the energy gaps and overlaps
frames = []
overlaps1 = []
overlaps2 = []
egaps1 = []
egaps2 = []
for i, orb in enumerate(orbs):
    fig, axes = plt.subplots(1, 3, figsize=(6.4*2, 4.8))
    # the first interaction
    Pd_5s = orb.sfos['Pd(5S)']
    sub_lumo = orb.sfos['Substrate(HOMO)']
    mos = [orb.mos['21AA'], orb.mos['17AA']]
    mix1 = pyfmo.analysis.mixing.Mixing(
        orb,
        sfos=[Pd_5s, sub_lumo], 
        mos=mos,
        connection_colors='b')

    overlaps1.append(abs(Pd_5s @ sub_lumo))
    egaps1.append(abs(Pd_5s.energy - sub_lumo.energy))


    # the second interaction
    Pd_4d = orb.sfos['Pd(2D:x2-y2)']
    sub_homo = orb.sfos['Substrate(LUMO)']
    mos = [orb.mos['22AA'], orb.mos['18AA']]
    mix2 = pyfmo.analysis.mixing.Mixing(
        orb,
        sfos=[Pd_4d, sub_homo], 
        mos=mos,
        connection_colors='r')

    overlaps2.append(abs(Pd_4d @ sub_homo))
    egaps2.append(abs(Pd_4d.energy - sub_homo.energy))


    # add the two mixings together
    main_mix = mix1 + mix2

    # draw the diagrams
    main_mix.draw_diagram(ax=axes[0], ylim=(-13, 2))
    plt.suptitle(f'C-H distance: {distances[i]: .2f} Å')

    # and then plot the overlaps
    axes[1].plot(distances[:i+1], overlaps1, label=r'$\langle 5S | 4AA \rangle$', c='b')
    axes[1].plot(distances[:i+1], overlaps2, label=r'$\langle 2D:x2-y2 | 5AA \rangle$', c='r')

    # and energy gaps
    axes[2].plot(distances[:i+1], egaps1, label=r'$| \varepsilon_{5S} - \varepsilon_{4AA} |$', c='b')
    axes[2].plot(distances[:i+1], egaps2, label=r'$| \varepsilon_{2D:x2-y2} - \varepsilon_{5AA} |$', c='r')

    # set some limits and labels
    axes[1].set_xlim(min(distances), max(distances))
    axes[1].set_ylim(0, 0.65)
    axes[1].vlines(1.614, 0, 0.6, colors='k', alpha=0.5, linestyle='dashed')
    axes[1].set_ylabel('Overlap')
    axes[1].set_xlabel('C-H distance / Å')

    axes[2].set_xlim(min(distances), max(distances))
    axes[2].set_ylim(0, 5)
    axes[2].vlines(1.614, 0, 5, colors='k', alpha=0.5, linestyle='dashed')
    axes[2].set_ylabel('Orbital Energy Gap / eV')
    axes[2].set_xlabel('C-H distance / Å')

    axes[1].legend(frameon=False)
    axes[2].legend(frameon=False)

    # save the frame and store the file location for later use
    frame_file = f'frames/{i}.png'
    frames.append(frame_file)

    plt.savefig(frame_file)
    plt.close()


# take the file locations of the frames we saved and make a movie
clips = [mvp.ImageClip(frame).set_duration(1/14) for frame in frames]
concat_clip = mvp.concatenate_videoclips(clips, method="compose")
concat_clip.write_videofile('orbint.mp4', fps=14)
concat_clip.write_gif('orbint.gif', fps=14)
