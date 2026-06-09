'''
This example shows how to find bonding and anti-bonding orbitals.

This example uses an ADF calculation on chloroethane
calculated at the OLYP/TZ2P level of theory.
The C-Cl bond was homolytically cleaved.
'''

import pyorbb


def bonding_score(mo, sfos1, sfos2):
	'''
	This function calculates the bonding score for an MO between two sets of SFOs.
	The bonding score is the sum of all products of coefficients and overlaps of the
	SFOs contributing to the MO.
	'''
	total = 0
	for sfo1 in sfos1:
		c1 = sfo1.coefficient(mo)
		for sfo2 in sfos2:
			c2 = sfo2.coefficient(mo)
			S = sfo1.overlap(sfo2)
			total += c1 * c2 * S

	return total

# load the orbital data
orbs = pyorbb.Orbitals('bonding.adf.rkf')

# and filter the correct orbitals
# for clarity we only show the filled alpha orbitals
mos = orbs.mos.filter(occupation=('fully_occupied', 'partially_occupied'), spin='A')

# split SFOs between the two fragments
sfos_fragment1 = orbs.sfos.filter(fragment='LeavingGroup')
sfos_fragment2 = orbs.sfos.filter(fragment='Substrate')

# go through each MO and calculate the score
rows = []
total_score = 0
for mo in mos:
	score = bonding_score(mo, sfos_fragment1, sfos_fragment2)
	total_score += score

	# if the score is too low we don't show it
	if abs(score) < 1e-3:
		continue

	# if the score is high enough we write a new string to be printed later
	rows.append([mo, score > 0, score])

# add a total row
rows.append(['Total', total_score > 0, total_score])

# print the gathered data
print('   MO  Bonding?  Score')
print('──────────────────────')
for mo, bonding, score in rows[:-1]:
	print(f'{str(mo):>5}  {str(bonding):8} {score: .3f}')
print('──────────────────────')
print(f'{str(rows[-1][0]):>5}  {str(rows[-1][1]):8} {rows[-1][2]: .3f}')
