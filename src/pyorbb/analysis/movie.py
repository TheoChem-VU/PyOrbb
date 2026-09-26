import pyorbb
import os
import numpy as np
import tcviewer
from typing import List
import moviepy as mvp
import uuid
import tcmu
from scipy.interpolate import RegularGridInterpolator


def make_frames(orbs: List[pyorbb.MO or pyorbb.FMO], 
        frames_dir: str = None, 
        transform: tcmu.geometry.Transform or List[tcmu.geometry.Transform] = None, 
        isovalue: float = 0.03, 
        color1=(1, 0, 0), 
        color2=(0, 0, 1)) -> List[str]:
    '''
    Generate frames for a movie showing the given orbitals and their molecules.

    Args:
        orbs: The orbital objects to draw.
        frames_dir: The directory to store the frames in.
        transform: The transformation to apply to the scenes containing the molecules and orbitals.
        isovalue: the isovalue for the isosurface drawn of the orbitals. We also draw the isosurface with the negative isovalue.
        color1: the color given to the positive isosurface.
        color2: the color given to the negative isosurface.

    Returns:
        A list of paths to the frames that were generated.
    '''

    # if no frames_dir is given we make a random one
    if frames_dir is None:
        # the dot at the start of the dir makes it hidden
        frames_dir = '.' + str(uuid.uuid4())

    # we generate the cube-files first
    cube_files = []
    for orb in orbs:
        cube_file = orb.cube_file()

        # orbital values often change signs across multiple frames
        # here we check if we should fix the sign or not
        if len(cube_files) > 0:
            # define a grid interpolator over the current orbital
            x, y, z = cube_file.x, cube_file.y, cube_file.z
            rgi = RegularGridInterpolator((x, y, z), cube_file.values.reshape(*cube_file.shape))

            # define the grid of the previous orbital
            x, y, z = cube_files[-1].x, cube_files[-1].y, cube_files[-1].z
            X, Y, Z = np.meshgrid(x, y, z)
            coords = np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).T
            # and compute the overlap between the current and previous orbital
            overlap = sum(cube_files[-1].values * rgi(coords))
            # if the sign was flipped we flip it back
            if overlap <= 0:
                cube_file.values *= -1

        cube_files.append(cube_file)
 
    # create a directory to store the frames in
    os.makedirs(frames_dir, exist_ok=True)

    frames = []
    # make a frame for every orbital object
    for i, cub in enumerate(cube_files):
        frame_file = os.path.join(frames_dir, f'{i}.png')
        frames.append(frame_file)

        # if we already have the frame somewhere we can skip
        if os.path.exists(frame_file):
            continue

        # open a new screen and draw the grid and molecule
        with tcviewer.Screen(headless=True) as scr:
            with scr.add_molscene() as scene:
                if isinstance(transform, list):
                    scene.transform = transform[i].to_vtkTransform()
                elif transform is not None:
                    scene.transform = transform.to_vtkTransform()
                scene.draw_molecule(cub.molecule)
                # we draw the positive and negative lobes separately
                scene.draw_isosurface(cub, isovalue=isovalue, color=color1, shininess=1, opacity=.5)
                scene.draw_isosurface(cub, isovalue=-isovalue, color=color2, shininess=1, opacity=.5)

            scene.screenshot(frame_file, enable_transparency=False)

    return frames


def make_orbital_movie(file: str, *args, fps: float = 60, **kwargs):
    '''
    Generate an mp4 movie showing the given orbitals and their molecules.

    .. seealso:: make_frames
    '''
    if os.path.exists(file):
        return

    if kwargs.get('frames_dir', None) is None:
        kwargs.setdefault('frames_dir', file +'.frames')

    frames = make_frames(*args, **kwargs)
    clips = [mvp.ImageClip(frame).with_duration(1/fps) for frame in frames]
    concat_clip = mvp.concatenate_videoclips(clips, method="compose")
    concat_clip.write_videofile(file, fps=fps)

