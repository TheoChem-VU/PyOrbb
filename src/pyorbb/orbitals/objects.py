"""Module defining the main classes used to access orbital data.

Data is organised hierarchically. The top-level |Orbitals| objects provide access to the |FMOs| and |MOs| objects, which provide access to individual |FMO| and |MO| objects.
The |Orbitals| objects serves as the loader of the calculation results.

Typical usage example:

  orbs = Orbitals('path/to/adf.rkf')
  

"""
import pyorbb
from scm import plams
import functools
import os
from typing import List, Dict, Tuple, Union
import math
import re

_ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else list(x)  # noqa: E731


class Orbital:
    '''
    Main class holding orbital information for |MO| and |FMO| objects.
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

        if self.singly_occupied or self.partially_occupied:
            if round(self.occupation) == 1:
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

        if self.singly_occupied or self.partially_occupied:
            if round(self.occupation) == 1:
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
        from tcmu.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, overwrite=overwrite, cube_file_prefix=cube_file_prefix) as job:
            [job.add_preamble(preamble) for preamble in preambles]
            job.rundir = os.path.split(self.parent.parent.kfpath)[0]
            job.name = 'densf'
            if isinstance(self, FMO):
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
        from tcmu.job.adf import DensfJob
        from tcintegral import grid
        import vtk

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, overwrite=overwrite) as job:
            [job.add_preamble(preamble) for preamble in preambles]
            job.generate_vtk()
            job.rundir = os.path.split(self.parent.parent.kfpath)[0]
            job.name = 'densf'
            if isinstance(self, FMO):
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

            if overwrite:
                os.remove(job.output_cub_paths[self])

            if job.can_skip():
                skipped = True
            else:
                skipped = False


        # add meta-data to the vtk file
        if not skipped:
            reader = vtk.vtkStructuredPointsReader()
            reader.SetFileName(job.output_cub_paths[self])
            reader.update()

            frag_mol = vtk.vtkFloatArray()
            frag_mol.SetName("Molecule")

            for atom in self.molecule:
                frag_mol.InsertNextValue(atom.atnum)
                frag_mol.InsertNextValue(atom.x)
                frag_mol.InsertNextValue(atom.y)
                frag_mol.InsertNextValue(atom.z)

            data = reader.GetOutput()
            data.GetFieldData().AddArray(frag_mol)

            writer = vtk.vtkStructuredPointsWriter()
            writer.SetFileName(job.output_cub_paths[self])
            writer.SetInputData(data)
            writer.Write()

        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_vtk_file(job.output_cub_paths[self])

    def draw(self, 
             gridsize: str = 'medium', 
             isovalue: float = 0.03, 
             overwrite: bool = False, 
             screen: "tcviewer.screen.Screen" = None,  # noqa: F821
             transform: "tcmu.geometry.Transform" = None):  # noqa: F821
        '''
        Generate and draw a cube-file for this |Orbital| object.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this |Orbital|.
            overwrite: whether to overwrite the previous calculation if found.
            screen: the ``tcviewer.screen.Screen`` object to use to draw this orbital. 
                If not given we start a new screen.
            transform: the geometrical tranfmormation to use with this orbital.

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
                   transform: "tcmu.geometry.Transform" = None) -> str:  # noqa: F821
        '''
        Generate a screenshot for this |Orbital| object.

        Args:
            output_path: the path to save the image to.
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this |Orbital|.
            overwrite: whether to overwrite the previous calculation if found.
            screen: the ``tcviewer.screen.Screen`` object to use to draw this orbital. 
                If not given we start a new screen.
            transform: the geometrical tranfmormation to use with this orbital.

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

    Each |MO| holds the following data that can be accessed like attributes.

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
        Calculate the total contribution of |FMO| objects from a specific fragment to this |MO|.
        The sum of all fragment characters is always ``1`` for each |MO|.

        Args:
            fragment: the fragment to calculate the character for.

        Example:

            .. code-block:: python

                >>> MO.fragment_character('NH3')
                0.469475215528633
                >>> MO.fragment_character('BH3')
                0.530524784471364
        '''
        fmos = self.parent.parent.fmos.filter(fragment=fragment)
        return sum(fmo.mulliken_contribution(self) for fmo in fmos)


    def coefficient(self, other: "FMO") -> float:
        '''
        Get the coefficient of an |FMO| into this |MO|.

        Args:
            other: the orbital that contributes to this |MO|.
        '''
        assert isinstance(other, FMO)

        return other.coefficient(self)


    def mulliken_contribution(self, other: "FMO") -> float:
        '''
        Get the Mulliken contribution of an |FMO| into this |MO|.

        Args:
            other: the orbital that contributes to this |MO|.
        '''
        assert isinstance(other, FMO)

        return other.mulliken_contribution(self)


class FMO(Orbital):
    '''
    Class holding data specifically for symmetry-adapted fragment orbitals.

    Each |FMO| holds the following data that can be accessed like attributes.

    .. list-table:: 
        :header-rows: 1

        * - Variable
          - Type
          - Description
        * - ``index``
          - ``int``
          - The index of this |FMO| in the overal |FMOs|.
        * - ``name``
          - ``str``
          - The regular name of this |FMO| as it would show up in ADFLevels.
        * - ``symmetry``
          - ``str``
          - The irreducible representation this |FMO| belongs to.
        * - ``symmetry_index``
          - ``int``
          - The index of this |FMO| in the overal |FMOs| that belong to the same irreducible representation.
        * - ``fragment``
          - ``str``
          - The name of the fragment the |FMO| belongs to.
        * - ``fragment_unique``
          - ``str``
          - If fragments do not have unique names (i.e. with atomic fragments) this name will be unique for the atom.
        * - ``fragment_index``
          - ``int``
          - The index of the |FMO| within the |FMOs| of the same fragment.
        * - ``spin``
          - ``str``
          - The spin of this |FMO|, either ``'A'``, ``'B'`` or ``'AB'``
        * - ``energy``
          - ``float``
          - The regular energy of the |FMO| in |kcal/mol|.
        * - ``approx_effective_energy``
          - ``float``
          - Approximated diagonal element of the Fock matrix belonging to the |FMO| in |kcal/mol|. This is available even if the Fock matrix cannot be read from the calculation.
        * - ``effective_energy``
          - ``float``
          - The diagonal element of the Fock matrix belonging to the |FMO| in |kcal/mol| if it could be read from the calculation.
        * - ``effective_energy_SCF0``
          - ``float``
          - The diagonal element of the Fock matrix after 0 SCF cycles belonging to the |FMO| in |kcal/mol| if it could be read from the calculation.
        * - ``occupation``
          - ``int``
          - The occupation number of this |FMO|. Either ``0``, ``1``, ``2``, or a fractional value if the electronic configuration is non-aufbau.
        * - ``occupied``
          - ``bool``
          - Whether the |FMO| has electrons in it.
        * - ``gross_population``
          - ``float``
          - The gross Mulliken population of this |FMO|.
        * - ``gross_spin``
          - ``float``
          - The gross Mulliken spin population of this |FMO|.
        * - ``molecule``
          - :class:`plams.Molecule`
          - The molecule object containing the atoms belonging to the fragment of this |FMO|.
    '''
    def __str__(self): 
        return self.make_name()

    def overlap(self, other: "FMO") -> float:
        '''
        Get the overlap between this |FMO| and another |FMO|.

        Args:
            other: the orbital to get the overlap with.

        .. note::

            The matmul operation ``@`` redirects to this method.
        '''
        assert isinstance(other, FMO)

        # these conditions apply due to orthonormality
        if self.spin != other.spin:
            return 0

        if self.symmetry != other.symmetry:
            return 0

        # access the right overlap matrix and return the right value
        S = self.parent.parent.data['matrices']['overlap'][self.symmetry][self.spin]
        return S[other.symmetry_index-1][self.symmetry_index-1]


    def fock(self, other: "FMO") -> float:
        '''
        Get the Fock matrix element between this |FMO| and another |FMO|.

        Args:
            other: the orbital to get the Fock matrix element with.
        '''
        assert isinstance(other, FMO)

        if self.spin != other.spin:
            return 0

        if self.symmetry != other.symmetry:
            return 0

        F = self.parent.parent.data['matrices']['fock'][self.symmetry][self.spin]
        return F[other.symmetry_index-1][self.symmetry_index-1]


    def mulliken_contribution(self, other: "MO", normalized=False) -> float:
        '''
        Get the mulliken contribution of this |FMO| into an |MO|.

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
        Get the coefficient of this |FMO| into an |MO|.

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


    def __matmul__(self, other: "FMO") -> float:
        '''
        Short-hand notation for getting the overlap with another |FMO|.
        '''
        return self.overlap(other)


    def make_name(self, spin: bool = True, frag_name: bool = True, relative_name: bool = False) -> str:
        '''
        Generate a name for this |FMO| with several options to modify it.

        Args:
            spin: whether to include spin in the name. It will be appended to the end as ``_{spin}``.
            frag_name: whether to include the fragment's unique name in the name as ``{fragment_unique}(...)``.
            relative_name: whether to use the relative name instead of the regular name.

        Examples:
            Generate the regular name of this |FMO|. This is the default name when printing the object.

            .. code-block:: python

                >>> fmo.make_name()
                'NH3(4A1)'
            
            One can also use relative naming.

            .. code-block:: python

                >>> fmo.make_name(relative_name=True)
                'NH3(LUMO)'

            One can also only get the name of the orbital by disabling the fragment name.

            .. code-block:: python

                >>> fmo.make_name(frag_name=False)
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
    Container class that stores information about both |MOs| and |FMOs|.
    |Orbitals| can also be given the paths to ``adf.rkf`` files from 
    related calculations to obtain more information. For example, the 
    path to a calculation with the number of SCF cycles set to 0 populates 
    the ``effective_energy_SCF0`` properties of the FMOs.

    Args:
        path: the path to an ``adf.rkf`` file containing information about the system of interest.
        path_SCF0: the path to an ``adf.rkf`` file containing information about a calculation with 0 SCF cycles.
            This argument is required to populate the ``FMO.effective_energy_scf0`` property
        path_fragments: dictionary containing fragment name as the key and path to its ``adf.rkf`` as the value.
        path_output: the path to an ``.out`` file generated by ADF. 
            This is required to read the kinetic energies for the MOs.

    Attributes:
        fmos (|FMOs|): the |FMOs| object storing the |FMO| objects associated with this system. Use this to select specific |FMO| for further analysis.
        mos (|MOs|): the |MOs| object storing the |MO| objects associated with this system.
        charges (Dict[str,int]): a dictionary storing formal charges of the complex and each fragment.

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
        self._gather_fmos()
        self._gather_mos()
        self._determine_formal_charges()
        self._gather_notices()


    def _get_data(self):
        self.data = pyorbb.orbitals.adf.read_data(self.reader, SCF0_reader=self.SCF0_reader, output=self.output)


    def _gather_fmos(self):
        self.fmos = FMOs([], self)
        fmo_mo_spin_match = self.data['calc_info']['unrestricted_mos'] == self.data['calc_info']['unrestricted_fmos']
        for fmo_idx in range(self.data['FMOs']['number']):
            for spin_idx, fmo_spin in enumerate(self.data['calc_info']['fmo_spins']):
                symlabel = self.data['FMOs']['symlabel'][fmo_idx]
                frag = self.data['FMOs']['fragment_types'][fmo_idx]

                if not fmo_mo_spin_match:
                    if self.data['calc_info']['unrestricted_mos']:
                        gross_pop = self.data['FMOs']['gross_population']['A'][fmo_idx] + self.data['FMOs']['gross_population']['B'][fmo_idx]
                        gross_spin = self.data['FMOs']['gross_population']['A'][fmo_idx] - self.data['FMOs']['gross_population']['B'][fmo_idx]
                    else:
                        gross_pop = self.data['FMOs']['gross_population']['AB'][fmo_idx]
                        gross_spin = 0
                else:
                    gross_pop = self.data['FMOs']['gross_population'][fmo_spin][fmo_idx]
                    gross_spin = 0

                data = {
                    'index': fmo_idx + 1 + self.data['MOs']['nfrozencores'][symlabel],
                    'name': self.data['FMOs']['adf_names'][fmo_spin][fmo_idx].removesuffix('_AB').removesuffix('_A').removesuffix('_B'),
                    'subspecies': self.data['FMOs']['subspecies'][fmo_idx],
                    'symmetry': symlabel,
                    'symmetry_index': self.data['FMOs']['symmetry_index'][fmo_idx] + 1 + self.data['MOs']['nfrozencores'][symlabel],
                    'densf_index': self.data['FMOs']['symmetry_index'][fmo_idx] + 1,
                    'fragment': frag,
                    # 'fragment_unique': self.data['FMOs']['fragment_unique']['total'][fmo_idx],
                    'fragment_index': self.data['FMOs']['fragment_index'][fmo_idx],
                    'spin': fmo_spin,
                    'energy': self.data['FMOs']['energy'][fmo_spin][fmo_idx] * 27.2114079527,
                    'occupation': float(self.data['FMOs']['occupation'][fmo_spin][fmo_idx]),
                    'occupied': int(self.data['FMOs']['occupation'][fmo_spin][fmo_idx]) > 0,
                    'gross_population': gross_pop,
                    'gross_spin': gross_spin,
                    'molecule': self.data['molecules'][frag],
                }
                if fmo_spin in self.data['FMOs']['approx_effective_energy']:
                    data['approx_effective_energy'] = self.data['FMOs']['approx_effective_energy'][fmo_spin][fmo_idx] * 27.2114079527
                else:
                    data['approx_effective_energy'] = (self.data['FMOs']['approx_effective_energy']['A'][fmo_idx]  + self.data['FMOs']['approx_effective_energy']['B'][fmo_idx]) * 27.2114079527

                if 'adf_names_fixed_principal' in self.data['FMOs']:
                    data['name'] = self.data['FMOs']['adf_names_fixed_principal'][fmo_spin][fmo_idx].removesuffix('_AB').removesuffix('_A').removesuffix('_B')

                # frag_unique = str(self.data['FMOs']['fragment_unique']['total'][fmo_idx])
                if symlabel in self.data['calc_info']['fmo_spinpolarizations'][frag]:
                    occs = self.data['calc_info']['fmo_spinpolarizations'][frag][symlabel]
                    spinpol = occs[0] - occs[1]
                    data['spin_pol'] = spinpol
                else:
                    data['spin_pol'] = 0

                data['effective_energy'] = None
                if 'effective_energy' in self.data['FMOs']:
                    data['effective_energy'] = self.data['FMOs']['effective_energy'][fmo_spin][fmo_idx] * 27.2114079527

                data['effective_energy_SCF0'] = None
                if 'effective_energy_SCF0' in self.data['FMOs']:
                    data['effective_energy_SCF0'] = self.data['FMOs']['effective_energy_SCF0'][fmo_spin][fmo_idx] * 27.2114079527

                fmo = FMO(data, self.fmos)
                self.fmos.orbitals.append(fmo)


    def _gather_mos(self):
        self.mos = MOs([], self)
        for moi in range(self.data['FMOs']['number']):
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
                fmo = MO(data, self.mos)
                self.mos.orbitals.append(fmo)


    def _gather_notices(self):
        self.notices = {'warning': [], 'error': [], 'info': []}
        if not self._check_effective_energies_available():
            self.notices['warning'].append(('No effective energies', 
'''Effective energies are not available for 
this calculation. To obtain them, please 
rerun the calculation with the following 
settings:

    Engine ADF
     PRINT FMATFMO
     FullFock Yes
     AllPoints Yes
    EndEngine
''', None))

        if any(c != 0 for c in self.charges.values()):
            self.notices['warning'].append(('Charged fragments', 
'''This system contains charged fragments.
We recommended you to check if effective
energies are required.
''', None))

        if self._check_spurious_mulliken_contr():
            fmos = self._mulliken_unstable_fmos()
            err = ''
            if len(fmos) > 0:
                max_fmo_len = max([len(str(fmo)) for fmo, _ in fmos])
                err += '\n  FMO'.ljust(max_fmo_len + 2) + '   Abs. Contr.\n'
                err += '  ' + '_' * (len('\n  FMO'.ljust(max_fmo_len + 2) + '   Abs. Contr.\n') - 4) + '\n'
                for fmo, score in fmos:
                    err += f'  {str(fmo):{max_fmo_len}} {score: .1%}\n'

            mos = self._mulliken_unstable_mos()
            if len(mos) > 0:
                max_mo_len = max([len(str(mo)) for mo, _ in mos])
                err += '\n  MO'.ljust(max_mo_len + 2) + '   Abs. Contr.\n'
                err += '  ' + '_' * (len('\n  MO'.ljust(max_mo_len + 2) + '   Abs. Contr.\n') - 4) + '\n'
                for mo, score in mos:
                    err += f'  {str(mo):{max_mo_len}} {score: .1%}\n'

            self.notices['warning'].append(('Mulliken artifacts', 
f'''We detected artifacts in the 
Mulliken analysis. Be carefull when 
interpreting Mulliken contributions, 
populations, and approximate effective 
energies from the following orbitals:
{err}''', [fmo for fmo, _ in fmos] + [mo for mo, _ in mos]))

        # check the EDA terms
        # OI: check irreps
        positive_eoi_irreps = []
        for lab in self.reader.read('Symmetry', 'symlab').split():
            l = lab.split(':')[0]
            if l not in positive_eoi_irreps:

                Eoi = float(self.reader.read('Energy', f'Orb.Int. {l}'))
                if Eoi > 0:
                    positive_eoi_irreps.append(l)

        if len(positive_eoi_irreps) > 0:
            s = r"    \n".join(positive_eoi_irreps)
            self.notices['error'].append(('Positive orb. int. energy', 
f'''The orbital interaction energy
is positive for the following irreps:
    {s}''', None))

        # check the electronic preparation
        for frag in self.fragments:
            fmos = [fmo for fmo in self.fmos if fmo.fragment == frag]
            polarized_fmos = self._polarized_fmos(fmos)
            if len(polarized_fmos) > 0:
                s = f'The "{frag}" fragment has at least\none large electronic shift\n\nMain polarized FMOs:\n'

                fmo_name_len = max([len(str(fmo)) for fmo, _ in polarized_fmos])
                for (fmo, dp) in polarized_fmos:
                    s += f'    {str(fmo):>{fmo_name_len}s}: {dp:+.2f} electrons\n'

                s += '\nCheck the electronic configuration!'

                self.notices['error'].append(('Incorrect electronic preparation', s, [r[0] for r in polarized_fmos]))

        if self._check_noninteger_occs():
            wrong_fmos = self._get_noninteger_occs()

            fmo_names = [str(fmo) for fmo in wrong_fmos]
            max_len = max(len(name) for name in fmo_names)
            occs = [f'{fmo.occupation:.2f}' for fmo in wrong_fmos]

            s = 'The following fractionally occupied\nFMOS were found:\n'
            for name, occ in zip(fmo_names, occs):
                s += f'    {name.ljust(max_len)} {occ} electrons\n'
            s += '\nCheck the electronic configuration!'

            self.notices['warning'].append(('Fractional occupations', s, wrong_fmos))


    def _polarized_fmos(self, fmos: List[FMO]) -> List[Tuple[FMO, float]]:
        '''
        Check given |FMO| objects for large changes in electronic population.

        Args:
            fmos: A list of |FMO| objects to check.

        Returns:
            A list of tuples with |FMO| objects and their difference in gross-population and occupation. The difference must be at least 0.7 electrons.
        '''
        polarized_fmos = []
        for fmo in fmos:
            dp = fmo.gross_population - fmo.occupation
            if abs(dp) > 0.7:
                polarized_fmos.append((fmo, dp))

        polarized_fmos = sorted(polarized_fmos, key=lambda r: -abs(r[1]))

        return polarized_fmos


    def _check_noninteger_occs(self):
        for fmo in self.fmos:
            if round(fmo.occupation) != fmo.occupation:
                return True

        return False


    def _get_noninteger_occs(self):
        ret = []
        for fmo in self.fmos:
            if round(fmo.occupation) != fmo.occupation:
                ret.append(fmo)

        return ret


    def _check_effective_energies_available(self):
        return any(hasattr(fmo, 'effective_energy') for fmo in self.fmos)


    def _mulliken_unstable_fmos(self):
        c = self.data['matrices']['mulliken_contribution']['total']
        ret = []
        for fmo in self.fmos:
            abs_C = sum(abs(c[:, fmo.index - 1]))
            if abs_C > 1.3:
                ret.append((fmo, abs_C))

        ret = sorted(ret, key=lambda row: -row[1])
        return ret


    def _mulliken_unstable_mos(self):
        c = self.data['matrices']['mulliken_contribution']['total']
        ret = []
        for mo in self.mos:
            abs_C = sum(abs(c[mo.index - 1]))
            if abs_C > 1.3:
                ret.append((mo, abs_C))

        ret = sorted(ret, key=lambda row: -row[1])
        return ret


    def _check_spurious_mulliken_contr(self):
        unstable_fmos = self._mulliken_unstable_fmos()
        unstable_mos = self._mulliken_unstable_mos()

        return len(unstable_fmos) > 0 or len(unstable_mos) > 0


    def _determine_formal_charges(self):
        # build up the effective charges of the atoms
        # this takes into account the atom number and number of frozen core electrons
        atomtypes = self.reader.read('Geometry', 'atomtype').split()
        eff_charges = self.reader.read('Geometry', 'atomtype effective charge')

        if isinstance(eff_charges, float):
            eff_charges = [eff_charges]

        if isinstance(atomtypes, float):
            atomtypes = [atomtypes]

        atomtype_charges = {typ: charge for typ, charge in zip(atomtypes, eff_charges)}

        # calculate the charges for the fragments and the complex
        charges = {}
        for frag in self.fragments:
            fmos = self.fmos.filter(fragment=frag)
            # we need the atoms in the molecule
            mol = fmos[0].molecule
            expected_Nelectrons = sum(atomtype_charges[atom.symbol] for atom in mol)
            actual_Nelectrons = round(sum(fmo.occupation for fmo in fmos))
            charges[frag] = expected_Nelectrons - actual_Nelectrons
            
        charges['complex'] = sum(charges.values())
        self.charges = charges


    @property
    def molecule(self) -> plams.Molecule:
        '''
        The molecule corresponding to the overall system.
        '''
        mol = plams.Molecule()
        for frag in self.fragments:
            fmo = [fmo for fmo in self.fmos if fmo.fragment == frag][0]
            mol += fmo.molecule
        return mol


    @property
    def fragments(self) -> List[str]:
        '''
        The names of the fragments defined in the calculation.
        '''
        return self.fmos.fragments


    def write_excel(self, out_file: str = None):
        '''
        Write the data corresponding to the system into an Excel file.

        Args:
            out_file: The filename of the Excel file to write.
        '''
        from pyorbb import write_excel

        if out_file is None:
            out_file = os.path.join(os.path.dirname(self.kfpath), 'PyOrbb.xlsx')
        write_excel.to_excel(self, out_file)


    @property
    def fmo_energy_types(self) -> List[str]:
        '''
        Get the orbital energy types that are available for the provided system.

        .. seealso:: 

            This property is a redirection of :attr:`FMOs.energy_types <pyorbb.orbitals.objects.FMOs.energy_types>`.
        '''
        return self.fmos.energy_types


    def rename_fragment(self, old: str, new: str):
        '''
        Rename the ``old`` fragment to ``new``.

        Args:
            old: the name of the fragment to rename.
            new: the name to rename the fragment to.

        Raises:
            ValueError: if the new name is already in use.

        .. seealso::

            See :attr:`Orbitals.fragments <pyorbb.orbitals.objects.Orbitals.fragments>` to obtain a list of fragment names that are currently used.
        '''
        if new in self.fragments:
            raise ValueError(f'Fragment name ``{new}`` is already in use.')

        for fmo in self.fmos:
            if fmo.fragment == old:
                fmo.fragment = new


    def get_mixer(self) -> "pyorbb.Mixer":
        return pyorbb.Mixer(self)



class OrbitalSelector:
    '''
    Class used to select |MOs| or |FMOs|. 
    It is responsible for decoding selection keys and filtering orbitals based on the selection key.

    Args:
        orbitals: a list of |FMOs| or |MOs| that will be managed by this class.
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

                >>> FMOs.get('NH3(HOMO)')
                NH3(3A1)
                >>> FMOs['NH3(HOMO)']
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

    def decode_key(self, key: Union[str, int]) -> dict:
        '''
        Decode a key into the relevant parts.
        Keys are given in the following format:

            {fragname}({orbname}[_{spin}][ {symmetry}])

        Where [:fragment_index], [_{spin}], and [ {symmetry}] are optional.

        If an FMO is desired you must begin the key with the fragment name
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

            Decode a key for an FMO specifying the fragment, orbname and spin.

            .. code-block:: python

                >>> FMOs.decode_key('NH3(1E1:1_B)')
                {'fragment': 'NH3', 'orbname': '1E1:1_B'}
    
            If multiple fragments have the same name (e.g. in a non-fragment analysis with atomic fragments)
            we can specify the fragment index with the colon.

            .. code-block:: python

                >>> FMOs.decode_key('C:4(1P:x)')
                {'fragment': 'C:4', 'orbname': '1P:x'}
        '''
        decoded = {
            'index': None,
            'global_index': None,
            'fragment': None,
            'orbname': None,
            'spin': None,
            'symmetry': None,
        }

        # if a single integer is given return only the index
        if isinstance(key, int):
            decoded['global_index'] = key
            return decoded

        # in case we have a FMO we need a fragment name
        fmo_regex = re.compile(r'(.+)\((\d+.+)\)_?([AB]?)')
        fmo_regex_result = fmo_regex.findall(key)
        if fmo_regex_result != []:
            decoded['fragment'], decoded['orbname'], decoded['spin'] = fmo_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = None
            return {k: v for k, v in decoded.items() if v is not None}

        # in case we have a FMO we need a fragment name
        fmo_relname_regex = re.compile(r'(.+)\(((?:HOMO|SOMO|LUMO|SUMO)(?:[+-]\d+)?)_?([AB])?\)')
        fmo_relname_regex_result = fmo_relname_regex.findall(key)
        if fmo_relname_regex_result != []:
            decoded['fragment'], decoded['orbname'], decoded['spin'] = fmo_relname_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = None
            return {k: v for k, v in decoded.items() if v is not None}

        # if the FMO regex fails we try the MO regex
        mo_regex = re.compile(r'(\d+[^_]+)_?([AB]?)')
        mo_regex_result = mo_regex.findall(key)
        if mo_regex_result != []:
            decoded['orbname'], decoded['spin'] = mo_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = None

            return {k: v for k, v in decoded.items() if v is not None}

        # if the FMO regex fails we try the MO regex
        mo_relname_regex = re.compile(r'((?:HOMO|SOMO|LUMO|SUMO)(?:[+-]\d+)?)_?([AB])?')
        mo_relname_regex_result = mo_relname_regex.findall(key)
        if mo_relname_regex_result != []:
            decoded['orbname'], decoded['spin'] = mo_relname_regex_result[0]
            if decoded['spin'] == '':
                decoded['spin'] = None

            return {k: v for k, v in decoded.items() if v is not None}

    def filter(self, 
            index: int or List[int] = None, 
            global_index: int or List[int] = None, 
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
            global_index: the global index of the orbital.
            symmetry: the symmetry label of the orbital.
            subspecies: the subspecies label of the orbital.
            spin: the spin label of the orbital, should be one of [``A``, ``B``, ``AB``].
            fragment: the fragment name of the FMO.
            fragment_index: the index of the fragment of the FMO.
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
            Select all FMOs of a given fragment.

            .. code-block:: python

                >>> FMOs.filter(fragment='NH3')
                [NH3(1A1), NH3(2A1), NH3(3A1), ...]

            Select all FMOs from the A2 irrep of the BH3 fragment.

            .. code-block:: python

                >>> FMOs.filter(symmetry='A2', fragment='BH3')
                [BH3(1A2), BH3(2A2), BH3(3A2), BH3(4A2)]

            Select all MOs that are named '1E1:1' or '1E1:2'.

            .. code-block:: python

                >>> MOs.filter(orbname=('1E1:1', '1E1:2'))
                [1E1:1, 1E1:2]

            Select the HOMO of the NH3 fragment.

            .. code-block:: python

                >>> FMOs.filter(orbname='HOMO', fragment='NH3')
                NH3(3A1)

            Get 1P orbitals for all carbons

            .. code-block:: python

                >>> FMOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C')
                [C:1(1P:x), C:1(1P:y), C:1(1P:z), C:2(1P:x), C:2(1P:y), C:2(1P:z), C:3(1P:x), C:3(1P:y), C:3(1P:z), C:4(1P:x), C:4(1P:y), C:4(1P:z)]
        
            Get 1P orbitals for the second carbon

            .. code-block:: python

                >>> FMOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C:2')
                [C:2(1P:x), C:2(1P:y), C:2(1P:z)]
                >>> FMOs.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C', fragment_index=2)
                [C:2(1P:x), C:2(1P:y), C:2(1P:z)]
        '''
        orbs = self.orbitals
        # filter down the orbitals in this object
        if index is not None:
            orbs = [orb for orb in orbs if orb.index in _ensure_list(index)]
        if global_index is not None:
            orbs = [orbs[i] for i in _ensure_list(global_index)]
        if symmetry is not None:
            orbs = [orb for orb in orbs if orb.symmetry in _ensure_list(symmetry)]
        if subspecies is not None:
            orbs = [orb for orb in orbs if orb.subspecies in _ensure_list(subspecies)]
        if spin is not None:
            orbs = [orb for orb in orbs if orb.spin in _ensure_list(spin)]

        # we match based on either fragment or fragment_unique
        # this ensures that if we select for instance "C(1P:x)" we match ALL carbons
        # if we match "C:1(1P:x)" we match only the first carbon
        if fragment is not None:
            orbs = [orb for orb in orbs if orb.fragment in _ensure_list(fragment)]

        if fragment_index is not None:
            orbs = [orb for orb in orbs if orb.fragment_index in _ensure_list(fragment_index)]

        # orbname can be either the proper name or the relative name
        if orbname is not None:
            orbs = [orb for orb in orbs if orb.name in _ensure_list(orbname) or orb.relative_name in _ensure_list(orbname)]

        # check for the occupation of the orbitals
        if occupation is not None:
            orbs_ = []
            for orb in orbs:
                for occ in _ensure_list(occupation):
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


class FMOs(OrbitalSelector):
    '''
    Object storing all |FMO| objects for the given calculation.
    '''
    @property
    def fragments(self) -> List[str]:
        '''
        Return a list of fragment names found in the orbitals.
        '''
        frags = []
        for fmo in self.orbitals:
            if fmo.fragment not in frags:
                frags.append(fmo.fragment)
        return frags

    @property
    def energy_types(self) -> List[str]:
        '''
        Object storing all |FMO| objects for the |Orbitals| objects.

        Returns:
            A list potentially containing ``energy``, ``effective_energy`` and ``effective_energy_SCF0``.
        '''
        ret = []
        if len(self.orbitals) > 0:
            orb = self.orbitals[0]
            if orb.energy is not None:
                ret.append('energy')
            if orb.effective_energy is not None:
                ret.append('effective_energy')
            if orb.approx_effective_energy is not None:
                ret.append('approx_effective_energy')
            if orb.effective_energy_SCF0 is not None:
                ret.append('effective_energy_SCF0')

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
    orbs = Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyOrbb/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.rkf')

    # for fmo in orbs.fmos:
    #     print(fmo.fragment_unique)

    # print(orbs.data.mos.kinetic_energy)

    # for mo in orbs.mos:
    #     print(mo, mo.kinetic_energy)
    # orbs.write_excel2()

    fmos = orbs.fmos.filter(symmetry='A2', fragment='BH3')
    print(fmos)
    fmos = orbs.fmos.filter(fragment='NH3')
    print(fmos)
    fmos = orbs.mos.filter(orbname=('1E1:1', '1E1:2'))
    print(fmos)
    fmos = orbs.fmos.filter(orbname='HOMO', fragment='NH3')
    print(fmos)

    mo = orbs.fmos.get('NH3(HOMO)')
    print(mo)

    mo = orbs.fmos['NH3(HOMO)']
    print(mo)

    mo = orbs.fmos['NH3']
    print(mo)

    dk = orbs.mos.decode_key('HOMO-2_A')
    print(dk)


    dk = orbs.fmos.decode_key('NH3(1E1:1_B)')
    print(dk)


    dk = orbs.fmos.decode_key('C:4(1P:x)')
    print(dk)

    orbs = Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/pyorbb/calculations/PyOrb_testing_2022/TransitionState/DielsAlder.Diene.results/adf.rkf')
    print(orbs.fmos.filter(orbname=('1P:x', '1P:y', '1P:z'), fragment='C', fragment_index=2))


    orbs = Orbitals('/Users/yumanhordijk/Downloads/pyr_c2v_frageda_occ2_pyridone.adf.rkf')
    for fmo in orbs.fmos.filter(fragment='CO'):
        if fmo.spin_pol != 0:
            print(fmo, fmo.spin_pol)
    for fmo in orbs.fmos.filter(fragment='NH'):
        if fmo.spin_pol != 0:
            print(fmo, fmo.spin_pol)
