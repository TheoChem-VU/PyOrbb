Warnings and Errors in PyOrbb
#############################

When loading a calculation into PyOrbb there may be a popup message mentioning a number of warnings and/or errors. This page provides more information about how these errors were detected and how they may potentially be remedied.

Errors
======

Errors are generally produced when the calculation setup is not correct.


Positive Orb. Int. Energy
-------------------------

Orbital interaction energies should always be a negative quantity (*i.e.* stabilizing). If it is positive it likely means that the electronic configuration was incorrectly set. Please check if the occupations are what you expect.


Incorrect Electronic Preparation
--------------------------------

PyOrbb will detect large electron gains or depletions to or from FMOs. PyOrbb will list FMOs with an absolute difference between the initial occupation and the Mulliken gross population larger than 0.7 electrons. If you receive this error check if the electron configurations of the fragments are correctly set.



Warnings
========

Warnings generally indicate recommendations for performing the analysis or issues with the system that are not easily fixed, but are moreso artifacts of the calculation.

Charged Fragments
-----------------

When PyOrbb detects that one or more of the fragments are charged it will warn you that the use of effective or approximated effective FMO energies may be important for your analysis. Generally, effective energies should yield more physically accurate analyses than using regular energies.

.. admonition:: Changing energy types
	:class: important

	To change the energy type of the FMOs, click the ``Energy Type`` button and select the energy type that you would like.

.. seealso::

	For more information, see the PyOrbb main article and Supporting Information Section S1.


No Effective Energies
---------------------

PyOrbb will warn you if it could not detect effective energies in the calculation. It will still provide you with the regular and approximate effective energies, which might be enough for your use-case.

.. admonition:: Calculating effective energies
	:class: important

	To obtain effective energies for your system; rerun the calculation using ADF and include the following in the input file::

		Engine ADF
		 PRINT FMATSFO
		 FullFock Yes
		 AllPoints Yes
		EndEngine

	PyOrbb should then have access to the effective energies.

.. seealso::

	For more information on the use of the energy types, see the PyOrbb main article and Supporting Information Section S1.


Mulliken Artifacts
------------------------

Mulliken analysis is known to be sensitive to systems with large orbital overlaps, and especially when large basis sets are employed. In those systems artifacts may appear in the form of negative Mulliken contributions. These negative contributions may lead to negative or larger than physically allowed Mulliken gross populations, and also affect the reliability of the approximate effective energies. PyOrbb identifies orbitals that are likely suffering from these artifacts by summing over the absolute contributions from or to the orbital.

Given the Mulliken contribution matrix (:math:`M \in \mathbb{R}^{N\times N}`) of the system and an FMO :math:`\psi_i` we give a warning if 

:math:`\sum_{l=1}^{N} |M_{il}| > 1.3`.

For an MO :math:`\Psi_k` we give a warning if

:math:`\sum_{j=1}^{N} |M_{jk}| > 1.3`.

The warning notice in the PyOrbb GUI will contain an overview of the affected orbitals. 


Fractional Occupations
----------------------

In some systems the FMOs may be fractionally occupied due to the symmetry of the fragment. For instance, this is very common when using atomic fragments. In reality, the FMOs should be occupied with an integer number of electrons. To remedy this, you may manually specify the occupations of the FMOs of the affected fragment.

.. seealso::
	
	`Click here <https://www.scm.com/doc/ADF/Input/Electronic_Configuration.html#orbital-occupations-electronic-configuration-excited-states>`_ for more information about occupation settings in ADF.
