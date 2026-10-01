
C–-C Bond Formation in Ethane
===============================

Usage
-----

The adf.rkf files for the uncatalysed and the Lewis-acid catalysed complexes scrutinized in this example, can be downloaded here:

1. **Download** :download:`Ethane.adf.rkf <./../../examples/ExamplesFromPaper/HomolyticEthane/Ethane.rkf>`
2. Open the PyOrbb application.
3. Drag the both .adf.rkf files to the GUI.

Now, the GUI will show two tabs. 
Each revealing the most important interactions. 
As the diene fragment is consistent in both complexes, one can select the orbital of interest (:math:`\pi`-HOMO\ :sub:`diene`\ --1) in both tabs to zoom in on what changes in the dienophile upon Lewis-acid coordination. 

4. Select the orbital of interest
5. Click the Filter button in the GUI



4. Select the orbital of interest.
5. Click on the orbital (:math:`\varphi`).

Now the PyOrbb viewer window will open with many options for drawing the orbitals.

Once the previous steps have been completed, the user should be able to reproduce the following analysis [1]_ [2]_.



Analysis
--------

The chemical bonding of molecules can be studied using Kohn-Sham MO analysis. 
In this section, we see that PyOrbb can assess the key interactions for the C--C bond formation and can readily isolate the orbital interactions within the pyramidalized H\ :sub:`3`\ C\ :sup:`•` fragment.

To investigate the C--C bond formation in ethane, the radical fragments were treated within the spin-unrestricted formalism. 
In
Figure 12a, the MO interaction scheme from the ADF-GUI is presented, for which the most important interactions were identified by PyOrbb and are summarized and schematically shown in the Figure below.
The key orbital interactions, including the pair bond interaction between the 3A\ :sub:`1` orbitals and the Pauli repulsive interaction between the 2A\ :sub:`1`, 1E\ :sub:`1`\ :sup:`(1)`, and 1E\ :sub:`1`\ :sup:`(2)` orbitals on each respective fragment, are correctly described as reported by Rodrigues Silva *et al.*


References
----------

.. [1] \ D. Rodrigues Silva, E. Blokker, J. M. Van Der Schuur, T. A. Hamlin, and F. M. Bickelhaupt, "Nature and strength of group-14 A–A′ bonds," Chemical Science 15 (2024): 1648–1656, https://doi.org/10.1039/d3sc06215e.
.. [2] |main art|