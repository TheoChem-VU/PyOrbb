import pyorbb
import matplotlib.pyplot as plt


orbs = pyorbb.Orbitals('../test/fixtures/NH3BH3/adf.rkf')
mixer = pyorbb.analysis.mixing.Mixer(orbs)
main_mix = pyorbb.analysis.mixing.Mixing(orbs)

for mix_ in mixer.orbital_interactions(N=4):
    main_mix += mix_

for mix_ in mixer.pauli_repulsions(N=2):
    main_mix += mix_

main_mix.sanitize()
main_mix.draw_diagram(simple=True)
plt.show()
