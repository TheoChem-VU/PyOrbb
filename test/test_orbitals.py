import pyorb
import os


def test_load_orbitals():
	orbs = pyorb.orbitals.Orbitals('fixtures/RadicalAddition/adf.rkf')


def test_energy1():
	orbs = pyorb.orbitals.Orbitals('fixtures/RadicalAddition/adf.rkf')
	assert round(orbs.mos['1A_A'].energy, 2) == -285.15


def test_energy2():
	orbs = pyorb.orbitals.Orbitals('fixtures/RadicalAddition/adf.rkf')
	assert round(orbs.mos['SOMO_A'].energy, 2) == -7.33
