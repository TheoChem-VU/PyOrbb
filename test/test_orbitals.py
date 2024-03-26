import pyorb
import os 

j = os.path.join


def test_load_orbitals():
	pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))


def test_energy1():
	orbs = pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['1A_A'].energy, 2) == -381.81


def test_energy2():
	orbs = pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['SOMO_A'].energy, 2) == -6.23


def test_energy3():
	orbs = pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO_A'] == orbs.mos['13A_A']


def test_energy4():
	orbs = pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO'] == orbs.mos['13A_A']


def test_energy5():
	orbs = pyorb.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SUMO'] == orbs.mos['13A_B']


if __name__ == '__main__':
	import pytest
	print(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	pytest.main()

