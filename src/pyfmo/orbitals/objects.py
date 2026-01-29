import pyfmo
from scm import plams
import functools
import os
from typing import List, Dict
import math
import platformdirs
import re

ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else list(x)  # noqa: E731


class Orbital:
    '''
    Main class holding orbital information for |MO| and |SFO| objects.
    This class is used to obtain information about the orbital, generate cube-files, and visualize orbitals.
    '''
    def __init__(self, data, parent):
        self.data = data
        for key, value in data.items():
            setattr(self, key, value)
        self.parent = parent

    def __repr__(self):
        return str(self)

    @functools.cached_property
    def relative_name(self) -> str:
        '''
        The relative name of the orbital. E.g. HOMO or HOMO-1
        '''
        orbitals = [orb for orb in self.parent.orbitals if orb.spin == self.spin and orb.spin_total_occupation == self.spin_total_occupation]
        if hasattr(self, 'fragment'):
            orbitals = [orb for orb in orbitals if orb.fragment == self.fragment]

        energies = sorted([orb.energy for orb in orbitals])
        order = energies.index(self.energy) + self.degeneracy_index

        if self.doubly_occupied:
            order = len(orbitals) - order - 1
            return f'HOMO-{order}' if order > 0 else 'HOMO'

        if self.singly_occupied:
            if self.occupation == 1:
                order = len(orbitals) - order - 1
                return f'SOMO-{order}' if order > 0 else 'SOMO'
            else:
                return f'SUMO+{order}' if order > 0 else 'SUMO'

        if self.unoccupied:
            return f'LUMO+{order}' if order > 0 else 'LUMO'

    @functools.cached_property
    def symmetry_relative_name(self) -> str:
        '''
        The relative name of the orbital in its irreducible representation. 
        E.g. the overall HOMO-2 could be the HOMO of its irreducible representation.
        '''
        orbitals = [orb for orb in self.parent.orbitals if orb.spin == self.spin and orb.spin_total_occupation == self.spin_total_occupation and orb.symmetry == self.symmetry]
        if hasattr(self, 'fragment'):
            orbitals = [orb for orb in orbitals if orb.fragment == self.fragment]

        energies = sorted([orb.energy for orb in orbitals])
        order = energies.index(self.energy) + self.degeneracy_index

        if self.doubly_occupied:
            order = len(orbitals) - order - 1
            return f'HOMO-{order}' if order > 0 else 'HOMO'

        if self.singly_occupied:
            if self.occupation == 1:
                order = len(orbitals) - order - 1
                return f'SOMO-{order}' if order > 0 else 'SOMO'
            else:
                return f'SUMO+{order}' if order > 0 else 'SUMO'

        if self.unoccupied:
            return f'LUMO+{order}' if order > 0 else 'LUMO'

    @functools.cached_property
    def fully_occupied(self) -> bool:
        '''
        Whether the orbital is fully occupied.
        '''
        if self.spin in ['A', 'B']:
            return self.singly_occupied
        if self.spin == 'AB':
            return self.doubly_occupied

    @functools.cached_property
    def partially_occupied(self) -> bool:
        '''
        Whether the orbital is not empty and not fully occupied.
        '''
        return not self.unoccupied and not self.fully_occupied

    @functools.cached_property
    def doubly_occupied(self) -> bool:
        '''
        Whether the orbital is doubly occupied.
        '''
        return self.spin_total_occupation == 2

    @functools.cached_property
    def singly_occupied(self) -> bool:
        '''
        Whether the orbital is singly occupied.
        '''
        return self.spin_total_occupation == 1

    @functools.cached_property
    def unoccupied(self) -> bool:
        '''
        Whether the orbital is unoccupied.
        '''
        return self.spin_total_occupation == 0

    @functools.cached_property
    def spin_total_occupation(self) -> int:
        '''
        The occupation of this orbital plus its spin counterpart if it exists.

        E.g. if orbital ``5A_A`` has an occupation of 1 and orbitals ``5A_B``has an 
        occupation of 0 then both orbitals will have the ``spin_total_occupation`` set to ``1``.
        '''
        return self.occupation + sum(orb.occupation for orb in self.spin_match_orbs)

    @functools.cached_property
    def spin_match_orbs(self) -> "Orbital":
        matching_orbs = [orb for orb in self.parent.orbitals if orb.name == self.name]
        if hasattr(self, 'fragment'):
            matching_orbs = [orb for orb in matching_orbs if orb.fragment == self.fragment]

        return [orb for orb in matching_orbs if orb != self]

    def cube_file(self, gridsize: str = 'medium', overwrite: bool = False, cube_file_prefix: str = None, preambles=[], grid_around_mol=None, gridextend=6):
        '''
        Generate a cube-file for this |Orbital| with a certain grid-size.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            overwrite: whether to overwrite the previous calculation if found.
            cube_file_prefix: prefix for the cube file path.

        .. seealso::
            :meth:`Orbital.draw` to draw and open a TCviewer screen showing this |Orbital|.
            :meth:`Orbital.screenshot` to generate a screenshot of this |Orbital|.
        '''
        from tcutility.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, overwrite=overwrite, cube_file_prefix=cube_file_prefix) as job:
            [job.add_preamble(preamble) for preamble in preambles]
            job.rundir = os.path.split(self.parent.parent.kfpath)[0]
            job.name = 'densf'
            if isinstance(self, SFO):
                job._sfos.append(self)
            else:
                job._mos.append(self)
            job.settings.ADFFile = self.parent.parent.kfpath
            job.cube_file_prefix = f"{self.parent.parent.kfpath}.densf/"
            os.makedirs(job.cube_file_prefix, exist_ok=True)
            if grid_around_mol is not None:
                job.grid_around_mol(grid_around_mol, extend=gridextend)
            else:
                job.gridsize(gridsize, extend=gridextend)


        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_cub_file(job.output_cub_paths[self])

    def vtk_file(self, gridsize: str = 'medium', overwrite: bool = False, preambles=[], grid_around_mol=None, gridextend=6):
        '''
        Generate a cube-file for this |Orbital| with a certain grid-size.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            overwrite: whether to overwrite the previous calculation if found.

        .. seealso::
            :meth:`Orbital.draw` to draw and open a TCviewer screen showing this |Orbital|.
            :meth:`Orbital.screenshot` to generate a screenshot of this |Orbital|.
        '''
        from tcutility.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, overwrite=overwrite) as job:
            [job.add_preamble(preamble) for preamble in preambles]
            job.generate_vtk()
            job.rundir = os.path.split(self.parent.parent.kfpath)[0]
            job.name = 'densf'
            if isinstance(self, SFO):
                job._sfos.append(self)
            else:
                job._mos.append(self)
            job.settings.ADFFile = self.parent.parent.kfpath
            job.cube_file_prefix = f"{self.parent.parent.kfpath}.densf/"
            os.makedirs(job.cube_file_prefix, exist_ok=True)
            if grid_around_mol is not None:
                job.grid_around_mol(grid_around_mol, extend=gridextend)
            else:
                job.gridsize(gridsize, extend=gridextend)

        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_vtk_file(job.output_cub_paths[self])

    def draw(self, 
             gridsize: str = 'medium', 
             isovalue: float = 0.03, 
             overwrite: bool = False, 
             screen: "tcviewer.screen.Screen" = None,  # noqa: F821
             transform: "tcutility.geometry.Transform" = None):  # noqa: F821
        '''
        Generate and draw a cube-file for this |Orbital| object.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this |Orbital|.
            overwrite: whether to overwrite the previous calculation if found.
            screen: the ``tcviewer.screen.Screen`` object to use to draw this orbital. 
                If not given we start a new screen.
            transform: the geometrical transformation to use with this orbital.

        .. seealso::
            :meth:`Orbital.cube_file` to generate and return a cube-file for this |Orbital|.
            :meth:`Orbital.screenshot` to generate a screenshot of this |Orbital|.
        '''
        import tcviewer

        # generate a cube-file or load an existing one
        cub = self.cube_file(gridsize=gridsize, overwrite=overwrite)

        if screen is None:
            scr = tcviewer.Screen()
            scr.__enter__()
            scr.window.show()
        else:
            scr = screen

        # and draw it with a specified isovalue
        with scr.add_molscene() as scene:
            c1, c2 = ([1, 0, 0], [0, 0, 1]) if self.occupied else ([1, .5, 0], [0, 1, 1])
            if transform is not None:
                scene.transform = transform.to_vtkTransform()

            scene.draw_molecule(self.molecule)
            scene.draw_isosurface(cub, -0.03, c1)
            scene.draw_isosurface(cub,  0.03, c2)

        if screen is None:
            scr.__exit__()

        return scr


    def screenshot(self, 
                   output_path: str = None, 
                   gridsize: str = 'medium', 
                   isovalue: float = 0.03, 
                   overwrite: bool = False, 
                   transform: "tcutility.geometry.Transform" = None) -> str:  # noqa: F821
        '''
        Generate a screenshot for this |Orbital| object.

        Args:
            output_path: the path to save the image to.
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this |Orbital|.
            overwrite: whether to overwrite the previous calculation if found.
            screen: the ``tcviewer.screen.Screen`` object to use to draw this orbital. 
                If not given we start a new screen.
            transform: the geometrical transformation to use with this orbital.

        .. seealso::
            :meth:`Orbital.cube_file` to generate and return a cube-file for this |Orbital|.
            :meth:`Orbital.draw` to draw and open a TCviewer screen showing this |Orbital|.
        '''
        import tcviewer

        if output_path is None:
            output_path = str(self) + '.png'

        # generate a cube-file or load an existing one
        cub = self.cube_file(gridsize=gridsize, overwrite=overwrite)

        # make a new screen and draw the orbital
        with tcviewer.Screen(headless=True) as scr:
            with scr.add_molscene() as scene:
                c1, c2 = ([1, 0, 0], [0, 0, 1]) if self.occupied else ([1, .5, 0], [0, 1, 1])
                if transform is not None:
                    scene.transform = transform.to_vtkTransform()

                scene.draw_molecule(self.molecule)
                scene.draw_isosurface(cub, -0.03, c1)
                scene.draw_isosurface(cub,  0.03, c2)

                # and take a screenshot
                scene.screenshot(output_path)

        return output_path


    @functools.cached_property
    def degeneracy_index(self) -> int:
        '''
        The index of this |Orbital| among its degenerate |Orbital| objects.
        '''
        return self.degenerate_orbitals.index(self)

    @functools.cached_property
    def degenerate(self) -> bool:
        '''
        Whether the |Orbital| is degenerate.
        '''
        return self.degeneracy > 1

    @functools.cached_property
    def degeneracy(self) -> int:
        '''
        The number of |Orbital| objects that are degenerate with this one.
        '''
        return len(self.degenerate_orbitals)

    @functools.cached_property
    def degenerate_orbitals(self) -> List["Orbital"]:
        '''
        |Orbital| objects that are very close in energy to this |Orbital|.
        '''
        return [orb for orb in self.parent.orbitals if math.isclose(orb.energy, self.energy, rel_tol=1e-8)]


class MO(Orbital):
    '''
    Class holding data specifically for molecular orbitals.

    Each |MO| holds the following data:

    .. list-table:: 
        :header-rows: 1

        * - Variable
          - Type
          - Description
        * - ``index``
          - ``int``
          - The index of this |MO| in the overal |MOs|.
        * - ``name``
          - ``str``
          - The regular name of this |MO| as it would show up in ADFLevels.
        * - ``symmetry``
          - ``str``
          - The irreducible representation this |MO| belongs to.
        * - ``symmetry_index``
          - ``int``
          - The index of this |MO| in the overal |MOs| that belong to the same irreducible representation.
        * - ``spin``
          - ``str``
          - The spin of this |MO|, either ``'A'``, ``'B'`` or ``'AB'``
        * - ``energy``
          - ``float``
          - The energy of the |MO| in |kcal/mol|.
        * - ``kinetic_energy``
          - ``float``
          - The kinetic energy of the |MO| in |kcal/mol| if it could be read from the calculation.
        * - ``occupation``
          - ``int``
          - The occupation number of this |MO|. Either ``0``, ``1`` or ``2``.
        * - ``occupied``
          - ``bool``
          - Whether the |MO| has electrons in it.
    '''
    def __str__(self):
        if self.spin == 'AB':
            return f'{self.name}'
        return f'{self.name}_{self.spin}'

    def fragment_character(self, fragment: str) -> float:
        '''
        Calculate the total contribution of |SFO| objects from a specific fragment to this |MO|.
        The sum of all fragment characters should be ``1`` for each |MO|.

        Args:
            fragment: the fragment to calculate the character for.

        Example:

            .. code-block:: python

                >>> MO.fragment_character('NH3')
                0.469475215528633
                >>> MO.fragment_character('BH3')
                0.530524784471364
        '''
        sfos = self.parent.parent.sfos.filter(fragment=fragment)
        return sum(sfo.mulliken_contribution(self) for sfo in sfos)


class SFO(Orbital):
    '''
    Class holding data specifically for symmetry-adapted fragment orbitals.

    Each |SFO| holds the following data:

    .. list-table:: 
        :header-rows: 1

        * - Variable
          - Type
          - Description
        * - ``index``
          - ``int``
          - The index of this |SFO| in the overal |SFOs|.
        * - ``name``
          - ``str``
          - The regular name of this |SFO| as it would show up in ADFLevels.
        * - ``symmetry``
          - ``str``
          - The irreducible representation this |SFO| belongs to.
        * - ``symmetry_index``
          - ``int``
          - The index of this |SFO| in the overal |SFOs| that belong to the same irreducible representation.
        * - ``fragment``
          - ``str``
          - The name of the fragment the |SFO| belongs to.
        * - ``fragment_unique``
          - ``str``
          - If fragments do not have unique names (i.e. with atomic fragments) this name will be unique for the atom.
        * - ``fragment_index``
          - ``int``
          - The index of the |SFO| within the |SFOs| of the same fragment.
        * - ``spin``
          - ``str``
          - The spin of this |SFO|, either ``'A'``, ``'B'`` or ``'AB'``
        * - ``energy``
          - ``float``
          - The energy of the |SFO| in |kcal/mol|.
        * - ``approx_site_energy``
          - ``float``
          - Approximated diagonal element of the Fock matrix belonging to the |SFO| in |kcal/mol|. This is available even if the Fock matrix cannot be read from the calculation.
        * - ``site_energy``
          - ``float``
          - The diagonal element of the Fock matrix belonging to the |SFO| in |kcal/mol| if it could be read from the calculation.
        * - ``site_energy_SCF0``
          - ``float``
          - The diagonal element of the Fock matrix after 0 SCF cycles belonging to the |SFO| in |kcal/mol| if it could be read from the calculation.
        * - ``occupation``
          - ``int``
          - The occupation number of this |SFO|. Either ``0``, ``1`` or ``2``.
        * - ``occupied``
          - ``bool``
          - Whether the |SFO| has electrons in it.
        * - ``gross_population``
          - ``float``
          - The gross Mulliken population of this |SFO|.
        * - ``gross_spin``
          - ``float``
          - The gross Mulliken spin population of this |SFO|.
        * - ``molecule``
          - :class:`plams.Molecule`
          - The molecule object containing the atoms belonging to the fragment of this |SFO|.
    '''
    def __str__(self): 
        return self.make_name()

    def overlap(self, other: "SFO") -> float:
        '''
        Get the overlap between this |SFO| and another |SFO|.

        Args:
            other: the orbital to get the overlap with.

        .. note::

            The matmul operation ``@`` redirects to this method.
        '''
        assert isinstance(other, SFO)

        # these conditions apply due to orthonormality
        if self.spin != other.spin:
            return 0

        if self.symmetry != other.symmetry:
            return 0

        # access the right overlap matrix and return the right value
        S = self.parent.parent.data['matrices']['overlap'][self.symmetry][self.spin]
        return S[other.symmetry_index-1][self.symmetry_index-1]


    def fock(self, other: "SFO") -> float:
        '''
        Get the Fock matrix element between this |SFO| and another |SFO|.

        Args:
            other: the orbital to get the Fock matrix element with.
        '''
        assert isinstance(other, SFO)

        if self.spin != other.spin:
            return 0

        if self.symmetry != other.symmetry:
            return 0

        F = self.parent.parent.data['matrices']['fock'][self.symmetry][self.spin]
        return F[other.symmetry_index-1][self.symmetry_index-1]


    def mulliken_contribution(self, other: "MO", normalized=False) -> float:
        '''
        Get the mulliken contribution of this |SFO| into an |MO|.

        Args:
            other: the orbital to get the Mulliken contribution with.
        '''
        assert isinstance(other, MO)

        if self.spin != other.spin and other.spin != 'AB':
            return 0

        if self.symmetry != other.symmetry:
            return 0

        if normalized:
            c = self.parent.parent.data['matrices']['mulliken_contribution_normalized'][self.symmetry][other.spin]
        else:
            c = self.parent.parent.data['matrices']['mulliken_contribution'][self.symmetry][other.spin]
        return c[other.symmetry_index-1][self.symmetry_index-1]


    def coefficient(self, other: "MO") -> float:
        '''
        Get the coefficient of this |SFO| into an |MO|.

        Args:
            other: the orbital to get the coefficient with.
        '''
        assert isinstance(other, MO)

        if self.spin != other.spin and other.spin != 'AB' and self.spin != 'AB':
            return 0

        if self.symmetry != other.symmetry:
            return 0

        c = self.parent.parent.data['matrices']['coefficients'][self.symmetry][other.spin]
        return c[other.symmetry_index-1][self.symmetry_index-1]


    def __matmul__(self, other: "SFO") -> float:
        '''
        Short-hand notation for getting the overlap with another |SFO|.
        '''
        return self.overlap(other)


    def make_name(self, spin: bool = True, frag_name: bool = True, relative_name: bool = False) -> str:
        '''
        Generate a name for this |SFO| with several options to modify it.

        Args:
            spin: whether to include spin in the name. It will be appended to the end as ``_{spin}``.
            frag_name: whether to include the fragment's unique name in the name as ``{fragment_unique}(...)``.
            relative_name: whether to use the relative name instead of the regular name.

        Examples:
            Generate the regular name of this |SFO|. This is the default name when printing the object.

            .. code-block:: python

                >>> sfo.make_name()
                'NH3(4A1)'
            
            One can also use relative naming.

            .. code-block:: python

                >>> sfo.make_name(relative_name=True)
                'NH3(LUMO)'

            One can also only get the name of the orbital by disabling the fragment name.

            .. code-block:: python

                >>> sfo.make_name(frag_name=False)
                '4A1'
        '''
        name = ''

        if frag_name:
            name += self.fragment + '('

        if relative_name:
            name += self.relative_name
        else:
            name += self.name

        if frag_name:
            name += ')'

        if self.spin != 'AB' and spin:
            name += f'_{self.spin}'

        return name

    @functools.cached_property
    def subspecies_relative_name(self) -> str:
        '''
        The relative name of the orbital in its irreducible representation. 
        E.g. the overall HOMO-2 could be the HOMO of its irreducible representation.
        '''
        orbitals = [orb for orb in self.parent.orbitals if orb.spin == self.spin and orb.spin_total_occupation == self.spin_total_occupation and orb.subspecies == self.subspecies]
        if hasattr(self, 'fragment'):
            orbitals = [orb for orb in orbitals if orb.fragment == self.fragment]

        energies = sorted([orb.energy for orb in orbitals])
        order = energies.index(self.energy) + self.degeneracy_index

        if self.doubly_occupied:
            order = len(orbitals) - order - 1
            return f'HOMO-{order}' if order > 0 else 'HOMO'

        if self.singly_occupied:
            if self.occupation == 1:
                order = len(orbitals) - order - 1
                return f'SOMO-{order}' if order > 0 else 'SOMO'
            else:
                return f'SUMO+{order}' if order > 0 else 'SUMO'

        if self.unoccupied:
            return f'LUMO+{order}' if order > 0 else 'LUMO'


class Orbitals:
    '''
    Container class that stores information about both |MOs| and |SFOs|.
    |Orbitals| can also be given the paths to ``adf.rkf`` files from 
    related calculations to obtain more information. For example, the 
    path to a calculation with the number of SCF cycles set to 0 populates 
    the ``site_energy_SCF0`` properties of the SFOs.

    Args:
        path: the path to an ``adf.rkf`` file containing information about the system of interest.
        path_SCF0: the path to an ``adf.rkf`` file containing information about a calculation with 0 SCF cycles.
            This argument is required to populate the ``SFO.site_energy_scf0`` property
        path_fragments: dictionary containing fragment name as the key and path to its ``adf.rkf`` as the value.
        path_output: the path to an ``.out`` file generated by ADF. 
            This is required to read the kinetic energies for the MOs.
    '''
    def __init__(self, path: str, path_SCF0: str = None, path_fragments: Dict[str, str] = None, path_output: str = None):
        self.reader = plams.KFReader(path)
        self.kfpath = os.path.abspath(path)
        self.SCF0_kfpath = path_SCF0
        self.SCF0_reader = plams.KFReader(path_SCF0) if path_SCF0 else None

        self.fragment_kfpaths = path_fragments
        if self.fragment_kfpaths:
            self.fragment_orbs = {frag: Orbitals(fpath) for frag, fpath in path_fragments.items()}
        else:
            self.fragment_orbs = {}

        self.output = os.path.abspath(path_output) if path_output else None

        self._get_data()
        self._gather_sfos()
        self._gather_mos()

    def _get_data(self):
        self.data = pyfmo.orbitals.adf.read_data(self.reader, SCF0_reader=self.SCF0_reader, output=self.output)

    def _gather_sfos(self):
        self.sfos = SFOs([], self)
        sfo_mo_spin_match = self.data['calc_info']['unrestricted_mos'] == self.data['calc_info']['unrestricted_sfos']
        for sfo_idx in range(self.data['SFOs']['number']):
            for spin_idx, sfo_spin in enumerate(self.data['calc_info']['sfo_spins']):
                symlabel = self.data['SFOs']['symlabel'][sfo_idx]
                frag = self.data['SFOs']['fragment_types'][sfo_idx]

                if not sfo_mo_spin_match:
                    if self.data['calc_info']['unrestricted_mos']:
                        gross_pop = self.data['SFOs']['gross_population']['A'][sfo_idx] + self.data['SFOs']['gross_population']['B'][sfo_idx]
                        gross_spin = self.data['SFOs']['gross_population']['A'][sfo_idx] - self.data['SFOs']['gross_population']['B'][sfo_idx]
                    else:
                        gross_pop = self.data['SFOs']['gross_population']['AB'][sfo_idx]
                        gross_spin = 0
                else:
                    gross_pop = self.data['SFOs']['gross_population'][sfo_spin][sfo_idx]
                    gross_spin = 0

                data = {
                    'index': sfo_idx + 1 + self.data['MOs']['nfrozencores'][symlabel],
                    'name': self.data['SFOs']['adf_names'][sfo_spin][sfo_idx].removesuffix('_AB').removesuffix('_A').removesuffix('_B'),
                    'subspecies': self.data['SFOs']['subspecies'][sfo_idx],
                    'symmetry': symlabel,
                    'symmetry_index': self.data['SFOs']['symmetry_index'][sfo_idx] + 1 + self.data['MOs']['nfrozencores'][symlabel],
                    'densf_index': self.data['SFOs']['symmetry_index'][sfo_idx] + 1,
                    'fragment': frag,
                    # 'fragment_unique': self.data['SFOs']['fragment_unique']['total'][sfo_idx],
                    'fragment_index': self.data['SFOs']['fragment_index'][sfo_idx],
                    'spin': sfo_spin,
                    'energy': self.data['SFOs']['energy'][sfo_spin][sfo_idx] * 27.2114079527,
                    'approx_site_energy': self.data['SFOs']['approx_site_energy'][sfo_spin][sfo_idx] * 27.2114079527,
                    'occupation': float(self.data['SFOs']['occupation'][sfo_spin][sfo_idx]),
                    'occupied': int(self.data['SFOs']['occupation'][sfo_spin][sfo_idx]) > 0,
                    'gross_population': gross_pop,
                    'gross_spin': gross_spin,
                    'molecule': self.data['molecules'][frag],
                }
                if 'adf_names_fixed_principal' in self.data['SFOs']:
                    data['name'] = self.data['SFOs']['adf_names_fixed_principal'][sfo_spin][sfo_idx].removesuffix('_AB').removesuffix('_A').removesuffix('_B')

                # frag_unique = str(self.data['SFOs']['fragment_unique']['total'][sfo_idx])
                if symlabel in self.data['calc_info']['sfo_spinpolarizations'][frag]:
                    occs = self.data['calc_info']['sfo_spinpolarizations'][frag][symlabel]
                    spinpol = occs[0] - occs[1]
                    data['spin_pol'] = spinpol
                else:
                    data['spin_pol'] = 0

                data['site_energy'] = None
                if 'site_energy' in self.data['SFOs']:
                    data['site_energy'] = self.data['SFOs']['site_energy'][sfo_spin][sfo_idx] * 27.2114079527

                data['site_energy_SCF0'] = None
                if 'site_energy_SCF0' in self.data['SFOs']:
                    data['site_energy_SCF0'] = self.data['SFOs']['site_energy_SCF0'][sfo_spin][sfo_idx] * 27.2114079527

                sfo = SFO(data, self.sfos)
                self.sfos.orbitals.append(sfo)

    def _gather_mos(self):
        self.mos = MOs([], self)
        for moi in range(self.data['SFOs']['number']):
            for spin_idx, mo_spin in enumerate(self.data['calc_info']['mo_spins']):
                symm_idx = self.data['MOs']['symmetry_index'][moi]
                symlabel = self.data['MOs']['symlabel'][moi]
                occ = int(self.data['MOs']['occupation'][symlabel][mo_spin][symm_idx])
                if 'kinetic_energy' in self.data['MOs']:
                    kin = self.data['MOs']['kinetic_energy'][symlabel][symm_idx] * 27.2114079527 if occ else 0
                else:
                    kin = None

                data = {
                    'index': moi + 1,
                    'name': f'{symm_idx+1}{symlabel}'.removesuffix('_AB').removesuffix('_A').removesuffix('_B'),
                    'symmetry': symlabel,
                    'symmetry_index': self.data['MOs']['symmetry_index'][moi] + 1,
                    'densf_index': self.data['MOs']['symmetry_index'][moi] + 1,
                    'spin': mo_spin,
                    'energy': self.data['MOs']['energy'][symlabel][mo_spin][symm_idx] * 27.2114079527,
                    'occupation': occ,
                    'occupied': int(self.data['MOs']['occupation'][symlabel][mo_spin][symm_idx]) > 0,
                    'kinetic_energy': kin,
                    'molecule': self.data['molecules']['complex'],
                }
                sfo = MO(data, self.mos)
                self.mos.orbitals.append(sfo)

    @property
    def molecule(self):
        mol = plams.Molecule()
        for frag in self.fragments:
            sfo = [sfo for sfo in self.sfos if sfo.fragment == frag][0]
            mol += sfo.molecule
        return mol

    @property
    def fragments(self):
        return self.sfos.fragments

    def write_excel(self, out_file: str = None):
        from pyfmo import write_excel

        if out_file is None:
            out_file = os.path.join(os.path.dirname(self.kfpath), 'pyfmo2.xlsx')
        write_excel.to_excel(self, out_file)

    @property
    def sfo_energy_types(self):
        return self.sfos.energy_types

    def rename_fragment(self, old: str, new: str):
        '''
        Rename the ``old`` fragment ``new``.
        '''
        for sfo in self.sfos:
            # we have to replace the ``fragment`` and ``fragment_unique`` properties
            if sfo.fragment == old:
                sfo.fragment = new
            # ``fragment_unique`` can contain ':'
            # if sfo.fragment_unique.split(':')[0] == old:
            #     if ':' not in sfo.fragment_unique:
            #         sfo.fragment_unique = new
            #         continue
            #     unique_part = sfo.fragment_unique.split(':')[1]
            #     sfo.fragment_unique = f'{new}:{unique_part}'


class OrbitalSelector:
    '''
    Class used to select |MOs| or |SFOs|. 
    It is responsible for decoding selection keys and filtering orbitals based on the selection key.

    Args:
        orbitals: a list of |SFOs| or |MOs| that will be managed by this class.
        parent: the parent |Orbitals| object.
    '''
    def __init__(self, orbitals: List[Orbital], parent: Orbitals):
        self.orbitals = orbitals
        self.parent = parent

    def __getitem__(self, key: int or str) -> List[Orbital] or Orbital:
        return self.get(key)

    def get(self, key: int or str) -> List[Orbital] or Orbital:
        '''
        Get |Orbital| objects based on the given key.

        Args:
            key: a string describing the orbital to be selected or the integer index of the orbital.

        Returns:
            A list of |Orbital| objects that match the given key.
            If there is only one return a single |Orbital| object.

        Examples:
            Select the HOMO of the NH3 fragment.

            .. code-block:: python

                >>> SFOs.get('NH3(HOMO)')
                NH3(3A1)
                >>> SFOs['NH3(HOMO)']
                NH3(3A1)

            Select a specific |MO|.

            .. code-block:: python

                >>> MOs.get('6A1')
                6A1
                >>> MOs['6A1']
                6A1

        .. seealso::

            :func:`~OrbitalSelector.filter` and :func:`~OrbitalSelector.decode_key`.

        .. note::

            The ``__getitem__`` method of this class redirects to this method, 
            allowing you to use indexing notation to obtain orbitals.
        '''
        return self.filter(**self.decode_key(key))

    def decode_key(self, key: str) -> dict:
        '''
        Decode a key into the relevant parts.
        Keys are given in the following format:

            {fragname}({orbname}[_{spin}][ {symmetry}])

        Where [:fragment_index], [_{spin}], and [ {symmetry}] are optional.

        If an SFO is desired you must begin the key with the fragment name
        and put the rest of the key within parentheses.

        Returns:
            A dictionary containing ``index``,  ``fragment``,
            ``orbname``, ``spin``, ``symmetry``.

        Examples:
            Decode a key specifying an MO.

            .. code-block:: python

                >>> MOs.decode_key('4A1')
                {'orbname': '4A1'}
            
            One can also use relative naming. Also specify alpha spin.

            .. code-block:: python

                >>> MOs.decode_key('HOMO-2_A')
                {'orbname': 'HOMO-2', 'spin': 'A'}

            Decode a key for an SFO specifying the fragment, orbname and spin.

            .. code-block:: python

                >>> SFOs.decode_key('NH3(1E1:1_B)')
                {'fragment': 'NH3', 'orbname': '1E1:1_B'}
    
            If multiple fragments have the same name (e.g. in a non-fragment analysis with atomic fragments)
            we can specify the fragment index with the colon.

            .. code-block:: python

                >>> SFOs.decode_key('C:4(1P:x)')
                {'fragment': 'C:4', 'orbname': '1P:x'}
        '''
        decoded = {
            'index': None,
            'fragment': None,
            'orbname': None,
            'spin': None,
            'symmetry': None,
        }

        # if a single integer is given return only the index
        if isinstance(key, int):
            decoded['index'] = key
            return decoded

        # in case we have a SFO we need a fragment name
        sfo_regex = re.compile(r'(.+)\((\d+.+)\)_?([AB]?)')
        sfo_regex_result = sfo_regex.findall(key)
        if sfo_regex_result != []:
            decoded['fragment'], decoded['orbname'], decoded['spin'] = sfo_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = 'AB'
            return {k: v for k, v in decoded.items() if v is not None}

        # in case we have a SFO we need a fragment name
        sfo_relname_regex = re.compile(r'(.+)\(((?:HOMO|SOMO|LUMO|SOMO)(?:[+-]\d+)?)(_[AB])?\)')
        sfo_relname_regex_result = sfo_relname_regex.findall(key)
        if sfo_relname_regex_result != []:
            decoded['fragment'], decoded['orbname'], decoded['spin'] = sfo_relname_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = 'AB'
            return {k: v for k, v in decoded.items() if v is not None}

        # if the SFO regex fails we try the MO regex
        mo_regex = re.compile(r'(\d+[^_]+)_?([AB]?)')
        mo_regex_result = mo_regex.findall(key)
        if mo_regex_result != []:
            decoded['orbname'], decoded['spin'] = mo_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = 'AB'

            return {k: v for k, v in decoded.items() if v is not None}

        # if the SFO regex fails we try the MO regex
        mo_relname_regex = re.compile(r'((?:HOMO|SOMO|LUMO|SOMO)(?:[+-]\d+)?)(_[AB])?')
        mo_relname_regex_result = mo_relname_regex.findall(key)
        if mo_relname_regex_result != []:
            decoded['orbname'], decoded['spin'] = mo_relname_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = 'AB'

            return {k: v for k, v in decoded.items() if v is not None}

    def filter(self, 
            index: int or List[int] = None, 
            symmetry: str or List[str] = None, 
            subspecies: str or List[str] = None, 
            spin: str or List[str] = None, 
            fragment: str or List[str] = None, 
            fragment_index: int or List[str] = None, 
            orbname: str or List[str] = None,
            occupation: float or str or List[float] or List[str] = None) -> Orbital or List[Orbital]:
        '''
        filter |Orbital| objects that match the given parameters.
        If any of the arguments is given as a ``Container`` we check for membership.

        Arguments:
            index: the index of the orbital.
            symmetry: the symmetry label of the orbital.
            subspecies: the subspecies label of the orbital.
            spin: the spin label of the orbital, should be one of [``A``, ``B``, ``AB``].
            fragment: the fragment name of the SFO.
            fragment_index: the index of the fragment of the SFO.
            orbname: the name of the orbital. Can be either the proper name or a relative name, 
                e.g. ``SOMO`` or ``LUMO+5``.
            occupation: what kind of occupation to allow. Can be a floating point number 
                specifying the occupation or a string from one of [``unoccupied``, ``partially_occupied``, ``fully_occupied``].
                Floating point numbers will be rounded to 2 decimals before comparison.

        Returns:
            The |Orbital| objects that match the provided arguments.
            If there is only one |Orbital| object selected, return only that one.
            Otherwise return a ``list`` of |Orbital| objects.
            Returns ``None`` if no matching |Orbital| objects were found.

        Examples:
            Select all SFOs of a given fragment.

            .. code-block:: python

                >>> SFOs.filter(fragment='NH3')
                [NH3(1A1), NH3(2A1), NH3(3A1), ...]

            Select all SFOs from the A2 irrep of the BH3 fragment.

            .. code-block:: python

                >>> SFOs.filter(symmetry='A2', fragment='BH3')
                [BH3(1A2), BH3(2A2), BH3(3A2), BH3(4A2)]

            Select all MOs that are named '1E1:1' or '1E1:2'.

            .. code-block:: python

                >>> MOs.filter(orbname=('1E1:1', '1E1:2'))
                [1E1:1, 1E1:2]

            Select the HOMO of the NH3 fragment.

            .. code-block:: python

                >>> SFOs.filter(orbname='HOMO', fragment='NH3')
                NH3(3A1)

            Get 1P orbitals for all carbons

            .. code-block:: python

                >>> SFOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C')
                [C:1(1P:x), C:1(1P:y), C:1(1P:z), C:2(1P:x), C:2(1P:y), C:2(1P:z), C:3(1P:x), C:3(1P:y), C:3(1P:z), C:4(1P:x), C:4(1P:y), C:4(1P:z)]
        
            Get 1P orbitals for the second carbon

            .. code-block:: python

                >>> SFOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C:2')
                [C:2(1P:x), C:2(1P:y), C:2(1P:z)]
                >>> SFOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C', fragment_index=2)
                [C:2(1P:x), C:2(1P:y), C:2(1P:z)]
        '''
        orbs = self.orbitals
        # filter down the orbitals in this object
        if index:
            orbs = [orb for orb in orbs if orb.index in ensure_list(index)]
        if symmetry is not None:
            orbs = [orb for orb in orbs if orb.symmetry in ensure_list(symmetry)]
        if subspecies is not None:
            orbs = [orb for orb in orbs if orb.subspecies in ensure_list(subspecies)]
        if spin is not None:
            orbs = [orb for orb in orbs if orb.spin in ensure_list(spin)]

        # we match based on either fragment or fragment_unique
        # this ensures that if we select for instance "C(1P:x)" we match ALL carbons
        # if we match "C:1(1P:x)" we match only the first carbon
        if fragment is not None:
            orbs = [orb for orb in orbs if orb.fragment in ensure_list(fragment)]

        if fragment_index is not None:
            orbs = [orb for orb in orbs if orb.fragment_index in ensure_list(fragment_index)]

        # orbname can be either the proper name or the relative name
        if orbname is not None:
            orbs = [orb for orb in orbs if orb.name in ensure_list(orbname) or orb.relative_name in ensure_list(orbname)]

        # check for the occupation of the orbitals
        if occupation is not None:
            orbs_ = []
            for orb in orbs:
                for occ in ensure_list(occupation):
                    if isinstance(occ, float):
                        if round(orb.occupation, 2) == round(occ, 2):
                            orbs_.append(orb)

                    if isinstance(occ, str):
                        if occ == 'occupied' and orb.occupied:
                            orbs_.append(orb)
                            continue
                        if occ == 'unoccupied' and orb.unoccupied:
                            orbs_.append(orb)
                            continue
                        if occ == 'fully_occupied' and orb.fully_occupied:
                            orbs_.append(orb)
                            continue
                        if occ == 'partially_occupied' and orb.partially_occupied:
                            orbs_.append(orb)
                            continue
            orbs = orbs_

        # return None if nothing was found
        if len(orbs) == 0:
            return None

        # squeeze the list if only one element exists
        if len(orbs) == 1:
            return orbs[0]

        return orbs

    def __len__(self):
        return len(self.orbitals)

    def __iter__(self):
        return iter(self.orbitals)

    @property
    def spins(self) -> List[str]:
        '''
        The spin species that are present in the given orbitals.
        '''
        return list(sorted({orb.spin for orb in self.orbitals}))

    @property
    def symmetry(self) -> List[str]:
        '''
        The spin species that are present in the given orbitals.
        '''
        return list(sorted({orb.symmetry for orb in self.orbitals}))

    @property
    def unrestricted(self) -> bool:
        '''
        Whether the calculation was performed in an unrestricted manner.
        '''
        return all(orb.spin in ['A', 'B'] for orb in self.orbitals)


class SFOs(OrbitalSelector):
    '''
    Object storing all |SFO| objects for the given calculation.
    '''
    @property
    def fragments(self) -> List[str]:
        '''
        Return a list of fragment names found in the orbitals.
        '''
        frags = []
        for sfo in self.orbitals:
            if sfo.fragment not in frags:
                frags.append(sfo.fragment)
        return frags

    @property
    def energy_types(self) -> List[str]:
        '''
        Object storing all |SFO| objects for the |Orbitals| objects.

        Returns:
            A list potentially containing ``energy``, ``site_energy`` and ``site_energy_SCF0``.
        '''
        ret = []
        if len(self.orbitals) > 0:
            orb = self.orbitals[0]
            if orb.energy is not None:
                ret.append('energy')
            if orb.site_energy is not None:
                ret.append('site_energy')
            if orb.approx_site_energy is not None:
                ret.append('approx_site_energy')
            if orb.site_energy_SCF0 is not None:
                ret.append('site_energy_SCF0')

        return ret

    @property
    def subspecies(self) -> List[str]:
        '''
        The spin species that are present in the given orbitals.
        '''
        return list(sorted({orb.subspecies for orb in self.orbitals}))


class MOs(OrbitalSelector):
    '''
    Object storing all |MO| objects for the |Orbitals| objects.
    '''
    ...



if __name__ == '__main__':
    orbs = Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.results/adf.rkf')

    # for sfo in orbs.sfos:
    #     print(sfo.fragment_unique)

    # print(orbs.data.mos.kinetic_energy)

    # for mo in orbs.mos:
    #     print(mo, mo.kinetic_energy)
    # orbs.write_excel2()

    sfos = orbs.sfos.filter(symmetry='A2', fragment='BH3')
    print(sfos)
    sfos = orbs.sfos.filter(fragment='NH3')
    print(sfos)
    sfos = orbs.mos.filter(orbname=('1E1:1', '1E1:2'))
    print(sfos)
    sfos = orbs.sfos.filter(orbname='HOMO', fragment='NH3')
    print(sfos)

    mo = orbs.sfos.get('NH3(HOMO)')
    print(mo)

    mo = orbs.sfos['NH3(HOMO)']
    print(mo)

    mo = orbs.sfos['NH3']
    print(mo)

    dk = orbs.mos.decode_key('HOMO-2_A')
    print(dk)


    dk = orbs.sfos.decode_key('NH3(1E1:1_B)')
    print(dk)


    dk = orbs.sfos.decode_key('C:4(1P:x)')
    print(dk)

    orbs = Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/TransitionState/DielsAlder.Diene.results/adf.rkf')
    print(orbs.sfos.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C', fragment_index=2))


    orbs = Orbitals('/Users/yumanhordijk/Downloads/pyr_c2v_frageda_occ2_pyridone.adf.rkf')
    for sfo in orbs.sfos.filter(fragment='CO'):
        if sfo.spin_pol != 0:
            print(sfo, sfo.spin_pol)
    for sfo in orbs.sfos.filter(fragment='NH'):
        if sfo.spin_pol != 0:
            print(sfo, sfo.spin_pol)
