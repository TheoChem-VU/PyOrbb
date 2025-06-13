import pyfmo
import itertools as it
import matplotlib.pyplot as plt


orbs = pyfmo.Orbitals('../test/fixtures/NH3BH3/adf.rkf')
mixer = pyfmo.analysis.mixing.Mixer(orbs)
main_mix = pyfmo.analysis.mixing.Mixing(orbs)

for mix_ in mixer.orbital_interactions(N=4):
    main_mix += mix_

for mix_ in mixer.pauli_repulsions(N=2):
    main_mix += mix_

main_mix.sanitize()
main_mix.draw_diagram(simple=True)
plt.show()
