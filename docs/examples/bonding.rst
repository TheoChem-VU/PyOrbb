
Bonding or Antibonding?
=======================

In this example we show how to differentiate between a bonding and antibonding MO for two chosen fragments.
As a system we chose the homolytic bond-cleavage of the carbon-chloride bond in chloroethane calculated at the ZORA-OLYP/TZ2P level of theory.

In general, to know if an MO has bonding or antibonding character for two chosen fragments we calculate the bond order

.. math::

   \textrm{BO}(\Psi_i) = \sum_{\psi_j^A \in \{\psi^A\}} \sum_{\psi_k^B \in \{\psi^B\}} \langle \psi_j^A | \Psi_i \rangle \langle \psi_k^B | \Psi_i \rangle \langle \psi_j^A | \psi_k^B \rangle,

where :math:`\Psi_i` is the MO of interest and :math:`\{\psi^A\}` and :math:`\{\psi^B\}` are the sets of SFOs from fragment A and B respectively. :math:`\langle \psi_j^A | \Psi_i \rangle` indicates the projection coefficient of :math:`\psi_j^A` onto the MO and :math:`\langle \psi_j^A | \psi_k^B \rangle` is the overlap between the SFOs. A negative value for :math:`\textrm{BO}(\Psi_i)` indicates antibonding character and a positive value bonding character. The magnitude of :math:`\textrm{BO}(\Psi_i)` indicates how strong the bond or antibond is. If the value is close to zero we are likely dealing with a non-bonding combination.

Calculating the bond order for the MOs of chloroethane yields the following bonding and antibonding MOs that have a bond order strength of 0.001 or greater:

.. tabs::

	.. tab:: Bonding

		.. list-table::
			:class: borderless

			* - .. image:: bonding_mos/resized_8A.png
			  - .. image:: bonding_mos/resized_11A.png
			  - .. image:: bonding_mos/resized_12A.png
			* - .. image:: bonding_mos/resized_13A.png
			  - .. image:: bonding_mos/resized_14A.png
			  - 

	.. tab:: Antibonding

		.. list-table::
			:class: borderless

			* - .. image:: bonding_mos/resized_9A.png
			  - .. image:: bonding_mos/resized_10A.png
			  - .. image:: bonding_mos/resized_16A.png
			* - .. image:: bonding_mos/resized_17A.png
			  -
			  - 


**Download** :download:`bonding.adf.rkf <../../examples/bonding_antibonding/bonding.adf.rkf>`

**Download** :download:`bonding.py <../../examples/bonding_antibonding/bonding.py>`


.. tabs:: 

	.. tab:: ``bonding.py`` 

		.. literalinclude:: ../../examples/bonding_antibonding/bonding.py
			:language: python

	.. tab:: Output

		The expected output is as follows

		.. code-block::

			MO      Bonding?   Score 
			──────────────────────────
			8A_A    True        0.091
			9A_A    False      -0.033
			10A_A   False      -0.004
			11A_A   True        0.026
			12A_A   True        0.036
			13A_A   True        0.028
			14A_A   True        0.014
			16A_A   False      -0.037
			17A_A   False      -0.044
			──────────────────────────
			Total   True        0.077

