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
import numpy as np
import moviepy.editor as mvp
import tcviewer
from typing import List
import cv2


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



def make_frames(orbs: List[pyfmo.Orbitals], sfo1_name: str, sfo2_name: str, frames_dir: str) -> List[str]:
    # create a directory to store the frames in
    os.makedirs(frames_dir, exist_ok=True)
    frames = []
    # make a frame for every orbital object
    for i, orb in enumerate(orbs):
        frame_file = os.path.join(frames_dir, f'{i}.png')
        frames.append(frame_file)

        # if we already have the frame somewhere we can skip
        if os.path.exists(frame_file):
            continue

        # obtain the SFOs
        sfo1 = orb.sfos[sfo1_name]
        sfo2 = orb.sfos[sfo2_name]

        # generate the cube-files
        # all rkf file are in the same directory,
        # we therefore pass a custom prefix to prevent
        # cube files from being overwritten
        cube_file_prefix = orb.kfpath.split('.')[0]
        cub1 = sfo1.cube_file(cube_file_prefix=cube_file_prefix)
        cub2 = sfo2.cube_file(cube_file_prefix=cube_file_prefix)

        # calculate the overlap cube file
        cub1.values *= cub2.values
        # to fix erroneous switching of the sign of the overlap grid
        # we multiply by the sign of the overlap given by ADF
        cub1.values *= (sfo1 @ sfo2) / abs(sfo1 @ sfo2)

        # open a new screen and draw the overlap grid and molecule
        with tcviewer.Screen(headless=True) as scr:
            with scr.add_molscene() as scene:
                scene.draw_molecule(orb.molecule)
                scene.draw_isosurface(cub1, isovalue=0.005, color=(1, 0, 1), shininess=1, opacity=.5)
                scene.draw_isosurface(cub1, isovalue=-0.005, color=(0, 1, 0), shininess=1, opacity=.5)

            # draw some informational text
            scene.draw_text(f'{sfo1_name}\n{sfo2_name}\nS = {abs(sfo1@sfo2): .2f}, Δϵ = {sfo1.energy - sfo2.energy: .1f} eV', fontsize=18)
            # generate a frame
            scene.screenshot(frame_file)

    return frames


# 4D/LUMO interaction
frames1 = make_frames(orbs, 'Pd(2D:x2-y2)', 'Substrate(5AA)', 'frames_int1')
# 5S/HOMO interaction
frames2 = make_frames(orbs, 'Pd(5S)', 'Substrate(4AA)', 'frames_int2')

# to make the movie we will concatenate the frames before writing a movie
os.makedirs('frames_final', exist_ok=True)
final_frames = []
for i, (frame1, frame2) in enumerate(zip(frames1, frames2)):
    # load the figures
    int1_fig = cv2.imread(frame1)
    int1_fig = cv2.resize(int1_fig, dsize=(480, 480), interpolation=cv2.INTER_CUBIC)
    int2_fig = cv2.imread(frame2)
    int2_fig = cv2.resize(int2_fig, dsize=(480, 480), interpolation=cv2.INTER_CUBIC)

    # and stitch them together
    final_fig = np.concatenate((int1_fig, int2_fig), axis=1)

    cv2.imwrite(f'frames_final/{i}.png', final_fig)
    final_frames.append(f'frames_final/{i}.png')


# take the file locations of the frames we saved and make a movie
clips = [mvp.ImageClip(frame).set_duration(1/14) for frame in final_frames]
concat_clip = mvp.concatenate_videoclips(clips, method="compose")
concat_clip.write_videofile('orbitals.mp4', fps=14)

clips = [mvp.ImageClip(frame).set_duration(1/14) for frame in frames1]
concat_clip = mvp.concatenate_videoclips(clips, method="compose")
concat_clip.write_gif('orbitals.gif', fps=14)

