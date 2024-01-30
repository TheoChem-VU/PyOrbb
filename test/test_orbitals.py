import pyfmo
import os 

j = os.path.join

## orbital load test -- radicals ##

def test_load_orbitals():
	pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))


def test_energy1():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['1A_A'].energy, 2) == -381.81


def test_energy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['SOMO_A'].energy, 2) == -6.23


def test_energy3():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO_A'] == orbs.mos['13A_A']


def test_energy4():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO'] == orbs.mos['13A_A']

def test_energy5():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SUMO'] == orbs.mos['13A_B']

## symmetry test for E' symmetries ##

def test_symm_orbitals():
	pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))

def test_symm_label():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.mos['HOMO-4'] == orbs.mos['EEE1:1']

def test_symm_energy1():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert round(orbs.mos['2 EEE1:2'].energy, 3) == -10.420


def test_symm_energy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert round(orbs.mos['HOMO-16'].energy, 3) == -14.603

if __name__ == '__main__':
	import pytest
	print(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	pytest.main()

