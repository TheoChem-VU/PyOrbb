
Basic Usage
===========

The main method to interact with orbital data in PyOrbb is with the |Orbitals| class.
This class allows you to load the data and also access it later on. 
Simply supply the ``adf.rkf`` file from your ADF calculation to the class to load the data.

.. code-block:: python
   
   >>> orbs = pyfmo.Orbitals('adf.rkf')

The data is then divided into the |MOs| and |SFOs| objects inside the |Orbitals| object.
The |MOs| and |SFOs| objects provide tools to select specific |MO| and |SFO| objects to analyse further.


Selecting Orbitals
------------------

To get an |MO| of the complex we can use the indexing notation of Python.

.. code-block:: python
   
   >>> homo = orbs.mos['2E1:2']
   >>> homo
   2E1:2

or using its relative name

.. code-block:: python

   >>> homo = orbs.mos['HOMO']
   >>> homo
   2E1:2

To select an |SFO| object we must also specify the fragment.

.. code-block:: python

   >>> nh3_sfo = orbs.sfos['NH3(LUMO)']
   >>> nh3_sfo
   NH3(4A1)

We can also select multiple objects at once using the :meth:`~pyfmo.orbitals.objects.OrbitalSelector.filter` methods of |MOs| and |SFOs|.
For example, to select all |MO| objects belonging to the A2 irreducible representation.

.. code-block::
   
   >>> a2_mos = orbs.mos.filter(symmetry='A2')
   >>> a2_mos
   [1A2, 2A2, 3A2, 4A2, 5A2, 6A2, 7A2, 8A2]

Or to select all |SFO| objects of the NH3 fragment.

.. code-block::
   
   >>> nh3_sfos = orbs.sfos.filter(fragment='NH3')
   >>> nh3_sfos
   [NH3(1A1), NH3(2A1), NH3(3A1), NH3(4A1), NH3(5A1), ...


Obtaining Data
--------------

The selected |MO| and |SFO| objects contain data pertaining to the orbitals themselves (e.g. orbital energies) but also to data relevant to the interaction of two orbitals (e.g. overlaps).
PyOrbb offers tools to easily obtain the required data.

To get the energy (in |eV|) of an |MO| or |SFO| we can access the ``energy`` attribute.

.. code-block::
   
   >>> homo.energy
   -6.392959708485001
   >>> nh3_sfo.energy
   -0.6265287931324357

Other properties can be accessed in the same way. See the |MO| and |SFO| documentation to see an overview of all properties available.

To obtain data related to the interaction between orbitals we need two objects. For example, to obtain the overlap between two |SFO| orbitals.

.. code-block::

   >>> nh3_homo = orbs.sfos['NH3(HOMO)']
   >>> bh3_lumo = orbs.sfos['BH3(LUMO)']
   >>> nh3_homo.overlap(bh3_lumo)
   -0.34332312250605734

Or to obtain the coefficient of an |SFO| into an |MO|.

.. code-block::

   >>> sfo = orbs.sfos['NH3(4A1)']
   >>> mo = orbs.mos['5A1']
   >>> sfo.coefficient(mo)
   0.02125011726149327


Advanced Examples
-----------------

For more advanced examples please see the `advanced examples <../examples/index.html>`_ section of this site.

.. toctree::

   ../examples/pi_orbitals
   ../examples/tracking
   ../examples/excitations
