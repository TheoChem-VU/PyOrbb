
Hetero- and Homolytic Bond Cleavage of Fluoromethane
=======================================================

Usage
-----

**Download** :download:`Homolytic.adf.rkf <./../../examples/ExamplesFromPaper/BondCleavage/Homo.rkf>`, :download:`Heterolytic_1.adf.rkf <./../../examples/ExamplesFromPaper/BondCleavage/Hetero1.rkf>`, and :download:`Heterolytic_2.adf.rkf <./../../examples/ExamplesFromPaper/BondCleavage/Hetero2.rkf>`




Analysis
--------

In the previous example, we investigated the homolytic bond cleavage of ethane and, subsequently, the fragmentation of the methyl radical. 
Interestingly, the same bond can be formed in more than one way, that is, from more than just the pair of fragments that is associated with the lowest-energy dissociation. 
This issue is addressed in the next example, which showcases the changes in orbital interaction mechanisms associated with the homolytic and heterolytic bond cleavage of fluoromethane using PyOrbb, following the work of Vermeeren and Bickelhaupt. 
The carbon-fluorine bond can be cleaved in three ways: homolytically, heterolytically with fluorine taking both bonding electrons, and heterolytically with CH\ :sub:`3` taking both bonding electrons. 
PyOrbb correctly identifies the correct orbital interaction mechanism for all three bond cleavage pathways.

In this example, we investigate charged species, which requires using the effective energies of the FMOs (Figure 14a, b and c) :sup:`9`. 
The effective energies refer to the FMO energies fully relaxed in the field of the molecule, as opposed to the orbital energies obtained from :math:`\boldsymbol{\psi}_{\mathrm{A}}` and :math:`\boldsymbol{\psi}_{\mathrm{B}}` of the isolated fragments (see Section S1). 
At this point, the FMO energies shift due to the inclusion of electrostatic potential between the fragments. 
We find a large shift in FMO energies when the effective energies are used for the two heterolytically cleaved systems (Figure 14b and e, and Figure 14c and f), indicating a strong stabilization of the FMOs of the anion due to the field of the cationic fragment, and in the same way, a strong destabilization of the FMOs of the cation. 
Conversely, the energies of the FMOs of the homolytically cleaved carbon-fluorine bond are negligibly affected by the field of the uncharged second fragment. 
The impact of energy choice is illustrated by the heterolytic cleavage of fluoromethane into CH\ :sub:`3`\ :sup:`+` and F\ :sup:`-` (Figure 14b and e). 
In Figure 14b, using effective orbital energies corrects the orbital interaction selection: the energy gap between the CH\ :sub:`3`\ :sup:`+`\ (4A') and F\ :sup:`-`\ (2A') FMOs (Figure 14e) is only 0.1 eV with regular orbital energies but 17.2 eV with effective energies, making the interaction appear much stronger under regular orbital energies.

References
----------

9. Author, A. B. Title of the paper. *Journal Name* **Year**,
   *Volume*, pages.
