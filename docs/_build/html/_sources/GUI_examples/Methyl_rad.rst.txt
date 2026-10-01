
Formation of the Methyl Radical
===============================

Usage
-----

The adf.rkf files for the uncatalysed and the Lewis-acid catalysed complexes scrutinized in this example, can be downloaded here:

1. **Download** :download:`Methyl_rad.adf.rkf <./../../examples/ExamplesFromPaper/HomolyticEthaneMethyl/Methyl.rkf>`
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

Once the previous steps have been completed, the user should be able to reproduce the following analysis [1]_ [2]_ .



Analysis
--------

To discuss the formation of the H\ :sub:`3`\ C\ :sup:`•` fragment, we fragment the pyramidal H\ :sub:`3`\ C\ :sup:`•` into the bare carbon atom, C\ :sup:`••••`, and the three hydrogen atoms, 3 H\ :sup:`•`. 
In Figure 13 we see that the pair bond formation occurs between the C(2s) and 1A\ :sub:`1`\ ', C(2p\ :sub:`x`) and the H\ :sub:`3`\ (1E\ :sub:`1`\ '\ :sup:`(1)`), and C(2p\ :sub:`y`) and the H\ :sub:`3`\ (1E\ :sub:`1`\ '\ :sup:`(2)`), respectively (for clarity, PyOrbb can separate the orbital mechanisms into the :math:`\alpha`- and :math:`\beta`-spin mechanisms, which are treated separately in the unrestricted formalism).
The C(2p\ :sub:`z`) is correctly excluded, as this is the H\ :sub:`3`\ C\ :sup:`•`\ (3A\ :sub:`1`) singly occupied orbital which partakes in pair bond formation in ethane, and, with a Mulliken contribution of 80.5%, is primarily formed by the C(2p\ :sub:`z`) orbital.



References
----------

.. [1] \ D. Rodrigues Silva, E. Blokker, J. M. Van Der Schuur, T. A. Hamlin, and F. M. Bickelhaupt, "Nature and strength of group-14 A–A′ bonds," Chemical Science 15 (2024): 1648–1656, https://doi.org/10.1039/d3sc06215e.
.. [2] |main art|



