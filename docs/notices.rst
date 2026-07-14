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


Mulliken Stability Score
------------------------

Mulliken analysis is known to be sensitive to systems with large orbital overlaps, and especially when large basis sets are employed. In those systems artifacts may appear leading to negative Mulliken gross populations and negative Mulliken contributions. PyOrbb identifies orbitals that are likely suffering from these artifacts by calculating a score for each 

PyOrbb calculates a score to measure the likelihood that these artifacts appeared in the provided calculation.

Given the Mulliken contribution matrix (:math:`C \in \mathbb{R}^{N\times N}`) of the system we calculate the stability score :math:`S` as:

:math:`S = \frac{1}{N} \sum_{i=1}^{N} \sum_{j=1}^{N} |C_{ij}|`.

By definition of the Mulliken analysis, the sum over the contribution matrix is exactly equal to :math:`N`. Therefore, in the ideal case that all contributions :math:`0 >= C_{ij} >= 1` the score :math:`S = 1`. Also due to the definition,any negative Mulliken contributions must be compensated by larger positive contributions and the sum over the absolute elements of the contribution matrix will be larger than :math:`N`. In that case the score :math:`S > 1` and we therefore detected artifacts in the Mulliken analysis. PyOrbb will give a warning when the score :math:`S >= 1.05`.

If PyOrbb gives a warning about the Mulliken analysis it will also provide the orbitals that are most affected. 


Fractional Occupations
----------------------

In some systems the FMOs may be fractionally occupied due to the symmetry of the fragment. For instance, this is very common when using atomic fragments. In reality, the FMOs should be occupied with an integer number of electrons. To remedy this, you may manually specify the occupations of the FMOs of the affected fragment.

.. seealso::
	
	`Click here <https://www.scm.com/doc/ADF/Input/Electronic_Configuration.html#orbital-occupations-electronic-configuration-excited-states>`_ for more information about occupation settings in ADF.
