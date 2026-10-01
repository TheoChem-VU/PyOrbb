
Diels-Alder Transition State
====================================================

Usage
-----

The adf.rkf file for this example can be downloaded here:

1. **Download** :download:`DielsAlder_TS.adf.rkf <./../../examples/ExamplesFromPaper/DielsAlder/DielsAlder.rkf>`
2. Open the PyOrbb application.
3. Drag the .adf.rkf file to the GUI.

Now the most important interactions will appear. To further study
fragment orbitals that interact, the user will want to visualize the
orbitals:

4. Select the orbital of interest.
5. Click on the orbital (:math:`\varphi`).

Now the PyOrbb viewer window will open with many options for drawing
the orbitals.

Once the previous steps have been completed, the user should be able
to reproduce the following analysis [1]_ [2]_.


Analysis
--------

Diels-Alder cycloadditions are a crucial class of reactions in
heterocyclic chemistry, materials chemistry, and total synthesis,
and have been extensively studied computationally. In this example,
we study the bonding mechanism between the diene, buta-1,3-diene,
and the dienophile, ethylene, in the Diels-Alder transition state.

Generally, in pericyclic transition states, two interactions play an
important role, namely the normal electron demand (NED) interaction
between a filled orbital on the diene and a virtual orbital of the
dienophile, and the inverse electron demand (IED), where electrons
are donated from a filled orbital on the dienophile to a virtual
orbital on the diene. PyOrbb correctly finds both the NED
and IED interactions (see below). Additionally, it identifies the
Pauli repulsive interaction between the HOMO-1 of the diene and the
HOMO of the dienophile, which contributes significantly to the Pauli
repulsion term in the EDA computation.


References
----------

.. [1] \ S. Yu, H. M. De Bruijn, D. Svatunek, T. A. Hamlin, and F. M. Bickelhaupt, "Factors Controlling the Diels–Alder Reactivity of Hetero‐1,3‐Butadienes," ChemistryOpen 7 (2018): 995–1004, https://doi.org/10.1002/open.201800193.
.. [2] |main art|