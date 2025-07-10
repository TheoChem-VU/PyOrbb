'''
This example shows how to analyse and filter excitation data from an ADF calculation.

The example system used is an electron donor-acceptor complex formed from indole and benzylbromide.
We are looking for transitions that have a wavelength higher than 300 nm, oscillator strength of
at least 1e-4 a.u.
The transitions corresponding to the excitations should have contributions higher than 10%
and the donating orbital for the transition should have at least 90% Donor fragment character
and the accepting orbital should have at least 90% Acceptor fragment character.
'''

import tcutility
import pyfmo

# exc should have wl of at least 300
wavelength_thresh          = 300
# and oscillator strength of at least 1e-4
oscillator_strength_thresh = 1e-4
# the contribution of a transition to an exc should be at least 10%
contribution_thresh        = 0.1
# the donor character of the donating MO should be at least 90%
donor_thresh               = 0.9
# and the acceptor character of the accepting MO should be at least 90%
acceptor_thresh            = 0.9

# we need to load the calculation with tcutility to read the excitation data
# it will be stored in the ``properties.excitations`` attribute.
res = tcutility.results2.read('excitations.adf.rkf')
# also load the orbital objects
orbs = pyfmo.Orbitals('excitations.adf.rkf')

# the excitation data is sorted by the irrep (only 'A' for this system) 
# and the excitation type (only 'SS' = singlet-singlet for this system)
excitation_data = res.properties.excitations
for irrep, irrep_data in excitation_data.items():
    for excitation_type, data  in irrep_data.items():
        for exc_index in range(data.number_of_excitations):
            wl = data.wavelengths[exc_index]
            # set a threshold for the wavelength
            if wl < wavelength_thresh:
                continue

            f12 = data.oscillator_strengths[exc_index]
            # and for oscillator strength
            if f12 < oscillator_strength_thresh:
                continue

            # obtain the contributions and MO names
            contribs = data.contributions[exc_index]
            from_MO_names = data.from_MO[exc_index]
            to_MO_names = data.to_MO[exc_index]

            n_transitions = len(contribs)
            transitions = []
            # each excitation has N transitions that contribute to it
            for transition_index in range(n_transitions):
                contrib = contribs[transition_index]
                # the contribution should be higher than 0.1
                if contrib < contribution_thresh:
                    continue

                # obtain the donating MO
                from_MO_name = from_MO_names[transition_index]
                from_MO = orbs.mos[from_MO_name]
                # the donor character should be high
                if from_MO.fragment_character('Donor') < donor_thresh:
                    continue

                # and the accepting MO
                to_MO_name = to_MO_names[transition_index]
                to_MO = orbs.mos[to_MO_name]
                # the acceptor character should be high
                if to_MO.fragment_character('Acceptor') < acceptor_thresh:
                    continue

                # if the transition survives we will write a string for it later
                transitions.append(f'{from_MO} ({from_MO.fragment_character("Donor"):.1%} Donor) -> {to_MO} ({to_MO.fragment_character("Acceptor"):.1%} Acceptor)')

            # write the data associated with this excitation
            print(f'Excitation[{irrep}, {excitation_type}, {exc_index+1}]:')
            print(f'  λ   = {wl:.1f} nm')
            print(f'  f12 = {f12:.4f} a.u.')
            print('  Transitions:')
            [print('    ' + transition) for transition in transitions]
            print()
