Oxidative Addition: C-H Bond Activation by Fe(CO)4
====================================================

Usage
-----

The adf.rkf file for this example can be downloaded here:

1. **Download** :download:`bonding.adf.rkf <../../examples/bonding_antibonding/bonding.adf.rkf>`
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
to reproduce the following analysis :sup:`1,2`.

Analysis
--------

Catalytic reactions are indispensable tools in modern synthetic
chemistry. The oxidative addition step is often the first and
rate-determining step in transition metal catalysis. In this
example, we analyze the transition state of the oxidative addition
of an iron complex into the methane C--H bond. PyOrbb correctly
identified the key orbital interaction mechanism (See below).

We see two major orbital interaction patterns, namely, the
interaction between CH\ :sub:`4`\ (5A) and Fe(CO)\ :sub:`4`\ (42A),
which corresponds to the :math:`\sigma`-donation from the C--H
orbital of CH\ :sub:`4` to the d orbital of Fe(CO)\ :sub:`4`. PyOrbb
also finds the :math:`\pi`-backdonation interaction between
CH\ :sub:`4`\ (6A) and Fe(CO)\ :sub:`4`\ (41A), which involves
donation of electrons from the d orbital of Fe(CO)\ :sub:`4` to the
:math:`\sigma^{*}`-C--H orbital of CH\ :sub:`4`.

References
----------

1. Author, A. B.; Author, C. D. Title of the paper. *Journal Name*
   **Year**, *Volume*, pages.
2. Author, E. F. Title of the paper or book. Publisher, City, Year.