'''
This example shows how to find bonding and anti-bonding orbitals.

This example uses an ADF calculation on chloroethane
calculated at the OLYP/TZ2P level of theory.
The C-Cl bond was homolytically cleaved.
'''

import pyorbb


def bond_order(mo, fmos1, fmos2):
	'''
	This function calculates the bond order for an MO between two sets of FMOS.
	The bonding order is the sum of all overlap populations between the two
	sets of FMOs.
	'''
	total = 0
	for fmo1 in fmos1:
		c1 = fmo1.coefficient(mo)
		for fmo2 in fmos2:
			c2 = fmo2.coefficient(mo)
			S = fmo1.overlap(fmo2)
			# this is the overlap population between fmo1 and fmo2
			total += 2 * mo.occupation * c1 * c2 * S

	return total

# load the orbital data
orbs = pyorbb.Orbitals('bonding.adf.rkf')

# split FMOS between the two fragments
fmos_fragment1 = orbs.fmos.filter(fragment='LeavingGroup')
fmos_fragment2 = orbs.fmos.filter(fragment='Substrate')

# to simplify we only consider one spin type for the MOs
mos = orbs.mos.filter(spin='A')

rows = []
total_bond_order = 0
# go through each MO and calculate the score
for mo in mos:
	score = bond_order(mo, fmos_fragment1, fmos_fragment2)
	total_bond_order += score

	# if the score is too low we don't show it
	if abs(score) < 1e-3:
		continue

	# if the score is high enough we write a new string to be printed later
	rows.append([mo, score > 0, score])

# add a total row
rows.append(['Total', total_bond_order > 0, total_bond_order])

# print the gathered data
print('   MO  Bonding?  Score')
print('──────────────────────')
for mo, bonding, score in rows[:-1]:
	print(f'{str(mo):>5}  {str(bonding):8} {score: .3f}')
print('──────────────────────')
print(f'{str(rows[-1][0]):>5}  {str(rows[-1][1]):8} {rows[-1][2]: .3f}')
