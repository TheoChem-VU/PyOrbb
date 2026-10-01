Lewis-Acid Catalysis: Diels-Alder
====================================================

Usage
-----

The adf.rkf files for the uncatalysed and the Lewis-acid catalysed complexes scrutinized in this example, can be downloaded here:

1. **Download** :download:`Uncatalysed_DA.adf.rkf <./../../examples/ExamplesFromPaper/LA_DielsAlder/TS_14.rkf>` and :download:`LA_catalysed_DA.adf.rkf <./../../examples/ExamplesFromPaper/LA_DielsAlder/TS_AlCl3_14.rkf>`
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

Once the previous steps have been completed, the user should be able to reproduce the following analysis [1]_ [2]_ [3]_.


Analysis
--------

The Diels-Alder reaction can be catalyzed by Lewis acids (LA).
Conventionally, LA catalysis has been explained as a consequence of enhanced reactivity resulting from stronger orbital interactions (*viz.* LUMO-lowering catalysis) :sup:`1b`. 
In 2020, Vermeeren *et al.* used ASM and EDA methodologies to point out that the enhanced reactivity upon LA-catalysis of the Diels-Alder reaction between isoprene and methyl acrylate was due to the reduced four-electron (Pauli) repulsion (*viz.* Pauli-lowering catalysis). 
They demonstrated that the LA pulls density away from the dienophile, resulting in a smaller amplitude of occupied :math:`\pi`-orbitals on the C--C double bond of the dienophile. 
This, in turn, results in less occupied-occupied overlap with the :math:`\pi`-system of the diene and thus less Pauli repulsion.

These insights can be obtained conveniently by applying the PyOrbb program to both the uncatalyzed and LA-catalyzed reaction (see below). 
To focus on the Pauli lowering effect, we hide the orbital interactions, leaving only the four-electron (Pauli) repulsive interactions in the interactive MO diagram. 
For both systems, the key interaction was identified to be between the occupied diene (:math:`\pi`-HOMO\ :sub:`diene`\ --1), which is the in-phase combination of the :math:`\pi`-system of the diene, and the :math:`\pi`-bonding orbital in the ester (:math:`\pi`-HOMO\ :sub:`dienophile`\  for both the uncatalyzed and the LA-catalyzed reaction). 
As reported in the literature, the lowering of the overlap between these two :math:`\pi`-systems upon LA coordination causes the reduced Pauli repulsion.


References
----------

.. [1] \ P. Vermeeren, T. A. Hamlin, I. Fernández, and F. M. Bickelhaupt, "How Lewis Acids Catalyze Diels–Alder Reactions," Angewandte Chemie 132 (2020): 6260–6265, https://doi.org/10.1002/ange.201914582; 
.. [2] \ T. A. Hamlin, F. M. Bickelhaupt, and I. Fernández, "The Pauli Repulsion-Lowering Concept in Catalysis," Accounts of Chemical Research 54 (2021): 1972–1981, https://doi.org/10.1021/acs.accounts.1c00016.
.. [3] |main art|
