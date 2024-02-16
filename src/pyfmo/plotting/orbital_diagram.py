import numpy as np
import matplotlib.pyplot as plt
import closed_interaction

# pip install matplotlib-label-lines
from labellines import labelLine, labelLines

import pyfmo

def degenarate_xvalues(Energies):
    
    # Setting the x-axis values for the scatters/energylevels     
    x = [0]
    for energy in range(len(Energies) - 1):
        current_energy = Energies[energy]
        next_energy = Energies[energy + 1]
        if int(current_energy) == int(next_energy):
            x[-1] = -0.2
            x.append(0.2)
        
        else:
            x.append(0)

    return x

def energy_yvalues(Orbitals):

    # puts the associated energies for MO or SFO in a list
    Energies = []
    for orbital in Orbitals:
        energy = orbital.energy
        Energies.append(energy)

    return Energies

def pair_occupation(Orbitals):

    # Puts the pair occupations for MO or SFO in a list
    Occupations = []
    for orbital in Orbitals:
        if orbital.occupation == 2:
                occupation = orbital.energy 
                Occupations.append(occupation)

    return Occupations

def single_occupation(Orbitals):
    
    # Puts the single occupations for MO or SFO in a list
    Occupations = []
    for orbital in Orbitals:
        if orbital.occupation == 1:
            if orb.spin == 'A':
                occupation = int(orbital.energy * 10) / 10
                Occupations.append(occupation)

    return Occupations

def occupation_xvalues(Occupations):

    x = [0]
    for energy in range(len(Occupations) - 1):
        current_energy = Occupations[energy]
        next_energy = Occupations[energy + 1]
        if int(current_energy) == int(next_energy):
            x[-1] = -0.2
            x.append(0.2)
        
        else:
            x.append(0)

    return x

def Energies_Contribution(MOs, SFOs):

    Energies = []  
    for MO in MOs:
        for SFO in SFOs:
            if orbitals.mulliken_contribution(SFO, MO) * 100 >= 1:
                SFO_energy = SFO.energy
                MO_energy = MO.energy
                Energies.append([MO_energy, SFO_energy])

    return Energies

def xvalues_Contributions(ECysF1, Cxy, CF1xy):

    result_list = []

    for pair in ECysF1:
        found = False
        for x1 in Cxy:
            if pair[0] == x1[0]:
                for x2 in CF1xy:
                    if pair[1] == x2[0]:
                        result_list.append([x1[1]+0.08, x2[1]-0.08])
                        found = True
                        break
                if found:
                    break

    return result_list

def percentages(MOs, SFOs):

    percentages = []
    for MO in MOs:
        for SFO in SFOs:
            percentage = orbitals.mulliken_contribution(SFO, MO) * 100
            if percentage >= 1:
                percentages.append(percentage)

    return percentages

def population(SFOs):

    populations = []
    for SFO in SFOs:
        population = sum([orbitals.mulliken_contribution(MO, SFO) * MO.occupation for MO in orbitals.mos])
        populations.append(population)

    return populations

def diagram(MOs, orbitals, SFOsF1, SFOsF2, color='b'):

    # Creating y-values
    # MO
    Energies_MO = energy_yvalues(MOs)
    Pair_Occupations_MO = pair_occupation(MOs)
    Single_Occupations_MO = single_occupation(MOs)

    # SFO1
    Energies_SFOs1 = energy_yvalues(SFOsF1)
    Pair_Occupations_SFOs1 = pair_occupation(SFOsF1)
    Single_Occupations_SFOs1 = single_occupation(SFOsF1)

    # SFO2
    Energies_SFOs2 = energy_yvalues(SFOsF2)
    Pair_Occupations_SFOs2 = pair_occupation(SFOsF2)
    Single_Occupations_SFOs2 = single_occupation(SFOsF2)


    # Creating x-values
    # MO
    x1 = degenarate_xvalues(Energies_MO)
    x2 = [(x + 1) for x in degenarate_xvalues(Energies_SFOs1)] 
    x3 = [(x - 1) for x in degenarate_xvalues(Energies_SFOs2)]

    # Pair SFO1 and SFO2
    x4 = occupation_xvalues(Pair_Occupations_MO)
    x5 = [(x + 1) for x in occupation_xvalues(Pair_Occupations_SFOs1)]
    x6 = [(x - 1) for x in occupation_xvalues(Pair_Occupations_SFOs2)]

    # Single SFO1 and SFO2
    x7 = occupation_xvalues(Single_Occupations_MO)
    x8 = [(x + 1) for x in occupation_xvalues(Single_Occupations_SFOs1)]
    x9 = [(x - 1) for x in occupation_xvalues(Single_Occupations_SFOs2)]


    # Plotting the scatters in the same plot
    # fig = plt.figure(figsize=(12,12))
    # ax1 = fig.add_subplot()
    ax1 = plt.gca()

    # MO
    ax1.scatter(x1, Energies_MO, s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x4, Pair_Occupations_MO, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if Single_Occupations_MO != []:
        x1.scatter(x7, Single_Occupations_MO, s=200, marker="$↑$", linewidth=0.2, c='k')

    # SFO1
    ax1.scatter(x2, Energies_SFOs1, s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x5, Pair_Occupations_SFOs1, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if Single_Occupations_SFOs1 != []:
        ax1.scatter(x8, Single_Occupations_SFOs1, s=200, marker="$↑$", linewidth=0.2, c='k')

    # SFO2
    ax1.scatter(x3, Energies_SFOs2,  s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x6, Pair_Occupations_SFOs2, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if Single_Occupations_SFOs2 != []:
        ax1.scatter(x9, Single_Occupations_SFOs2, s=200, marker="$↑$", linewidth=0.5, c='k')

    # Plotting the lines
    Cxy = [list(x) for x in zip(Energies_MO, x1)]
    CF1xy = [list(x) for x in zip(Energies_SFOs1, x2)]
    CF2xy = [list(x) for x in zip(Energies_SFOs2, x3)]

    # generate the x and y-values
    ECysF1 = Energies_Contribution(MOs, SFOsF1)
    ECysF2 = Energies_Contribution(SFOsF2, MOs)
    ECxsF1 = xvalues_Contributions(ECysF1, Cxy, CF1xy)
    ECxsF2 = xvalues_Contributions(ECysF2, CF2xy, Cxy)

    # Generate percentages to label the lines
    percentageF1 = percentages(MOs, SFOsF1)
    percentageF2 = percentages(MOs, SFOsF2)

    # Plotting the lines
    xvals = []
    for xpair, ypair, percentage in zip(ECxsF1, ECysF1, percentageF1):
        ax1.plot(xpair, ypair, alpha=0.5, ls='dashed', c=color, label=f"{int(percentage)/10 * 10} %")
        xvals.append(0.5)

    for xpair2, ypair2, percentage2 in zip(ECxsF2, ECysF2, percentageF2):
        ax1.plot(xpair2, ypair2, alpha=0.5, ls='dashed', c=color, label=f"{int(percentage2)/10 * 10} %")
        xvals.append(-0.5)

    # inline label for the contribution
    labelLines(ax1.get_lines(), align=True, xvals=xvals, ha="center", va="center", backgroundcolor="none", fontsize=6)


    # Getting the populations
    popF1 = population(SFOsF1)
    popF2 = population(SFOsF2)


    # Plotting the orbital label and energies for MO, SFO1 and SFO2
    # MO
    for i, j, l in zip(x1, Energies_MO, MOs):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(40,4), size=7, ha="center", va="top", textcoords="offset points")

    # SFO1
    for i, j, l, pop in zip(x2, Energies_SFOs1, SFOsF1, popF1):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(40,8), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")
        ax1.annotate(f"{round(pop, 2)} e", xy=(i, j), xytext=(40,-2), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")

    # SFO2
    for i, j, l, pop in zip(x3, Energies_SFOs2, SFOsF2, popF2):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(-40,8), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")
        ax1.annotate(f"{round(pop, 2)} e", xy=(i, j), xytext=(-40,-2), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")

    #Tuning Axis
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.xticks([-1, 0, 1], ['SFO1', "MO", 'SFO2'])
    plt.ylabel(r'$\epsilon$ / eV')
    
    # Plot
    plt.xlim(-1.7, 1.7)
    # plt.show()



if __name__ == '__main__':
    fig = plt.figure(figsize=(12,12))
    orbitals = pyfmo.orbitals.Orbitals('/Users/Tori/PyFMO/test/fixtures/NH3BH3/adf.rkf')
    MOs = orbitals.mos['HOMO-3':'LUMO+3']
    SFOsF1 = orbitals.sfos['Donor(HOMO-2)':'Donor(LUMO+2)']
    SFOsF2 = orbitals.sfos['Acceptor(HOMO-2)':'Acceptor(LUMO+2)']
    mixing2 = max(closed_interaction.get_2mixings(orbitals, orbitals.mos['HOMO-2']))
    print(mixing2)
    fragments = list(orbitals.fragments)
    diagram(MOs, orbitals, SFOsF1, SFOsF2)
    diagram(mixing2.mos, orbitals, [sfo for sfo in mixing2.sfos if sfo.fragment == fragments[0]], [sfo for sfo in mixing2.sfos if sfo.fragment == fragments[1]], color='r')
    plt.show()

