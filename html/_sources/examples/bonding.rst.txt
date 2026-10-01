
Bonding or Antibonding?
=======================

In this example we show how to differentiate between a bonding and antibonding MO for two chosen fragments.
As a system we chose the homolytic bond-cleavage of the carbon-chloride bond in chloroethane calculated at the ZORA-OLYP/TZ2P level of theory.

In general, to know if an MO has bonding or antibonding character for two chosen fragments we calculate the bond order

.. math::

   \textrm{BO}(\Psi_i) = 2N_i \sum_{\psi_j^A \in \{\psi^A\}} \sum_{\psi_k^B \in \{\psi^B\}} \langle \psi_j^A | \Psi_i \rangle \langle \psi_k^B | \Psi_i \rangle \langle \psi_j^A | \psi_k^B \rangle,

where :math:`\Psi_i` is the MO of interest with occupation :math:`N_i` and :math:`\{\psi^A\}` and :math:`\{\psi^B\}` are the sets of SFOs from fragment A and B respectively. :math:`\langle \psi_j^A | \Psi_i \rangle` is the projection coefficient of :math:`\psi_j^A` onto the MO and :math:`\langle \psi_j^A | \psi_k^B \rangle` is the overlap between the SFOs. A negative value for :math:`\textrm{BO}(\Psi_i)` is antibonding character and a positive value bonding character. The magnitude of :math:`\textrm{BO}(\Psi_i)` indicates how strong the bond or antibond is. If the value is close to zero we are likely dealing with a non-bonding combination.

Calculating the bond order for the MOs of chloroethane yields the following bonding and antibonding MOs that have a bond order strength of 0.001 or greater:

.. tabs::

	.. tab:: Bonding

		.. list-table::
			:class: borderless

			* - .. image:: bonding_mos/4A.png
			  - .. image:: bonding_mos/resized_8A.png
			  - .. image:: bonding_mos/resized_11A.png
			* - MO 4A
			  - MO 8A
			  - MO 11A
			* - .. image:: bonding_mos/resized_12A.png
			  - .. image:: bonding_mos/resized_13A.png
			  - .. image:: bonding_mos/resized_14A.png
			* - MO 12A
			  - MO 13A
			  - MO 14A

	.. tab:: Antibonding

		.. list-table::
			:class: borderless

			* - .. image:: bonding_mos/resized_9A.png
			  - .. image:: bonding_mos/resized_10A.png
			  - .. image:: bonding_mos/resized_16A.png
			* - MO 9A
			  - MO 10A
			  - MO 16A
			* - .. image:: bonding_mos/resized_17A.png
			  -
			  - 
			* - MO 17A
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

			   MO  Bonding?  Score
			──────────────────────
			 4A_A  True      0.001
			 8A_A  True      0.183
			 9A_A  False    -0.067
			10A_A  False    -0.009
			11A_A  True      0.052
			12A_A  True      0.072
			13A_A  True      0.056
			14A_A  True      0.027
			16A_A  False    -0.073
			17A_A  False    -0.088
			──────────────────────
			Total  True      0.154

