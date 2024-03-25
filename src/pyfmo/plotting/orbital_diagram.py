import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines
import pyfmo

def degenerate_xvalues(energies):
    
    # Setting the x-axis values for the scatters/energylevels     
    x = [0]
    for energy in range(len(energies) - 1):
        current_energy = energies[energy]
        next_energy = energies[energy + 1]
        if int(current_energy) == int(next_energy):
            x[-1] = -0.2
            x.append(0.2)
        
        else:
            x.append(0)

    return x

def energy_yvalues(orbitals):

    # puts the associated energies for MO or SFO in a list
    energies = []
    for orbital in orbitals:
        energy = orbital.energy
        energies.append(energy)

    return energies

def pair_occupation(orbitals):

    # Puts the pair occupations for MO or SFO in a list
    occupations = []
    for orbital in orbitals:
        if orbital.occupation == 2:
                occupation = orbital.energy 
                occupations.append(occupation)

    return occupations

def single_occupation(orbitals):
    
    # Puts the single occupations for MO or SFO in a list
    occupations = []
    for orbital in orbitals:
        if orbital.occupation == 1:
            if orbital.spin == 'A':
                occupation = int(orbital.energy * 10) / 10
                occupations.append(occupation)

    return occupations

def occupation_xvalues(occupations):

    x = [0]
    for energy in range(len(occupations) - 1):
        current_energy = occupations[energy]
        next_energy = occupations[energy + 1]
        if int(current_energy) == int(next_energy):
            x[-1] = -0.2
            x.append(0.2)
        
        else:
            x.append(0)

    return x

def energies_contribution(mos, sfos):

    energies = []  
    for mo in mos:
        for sfo in sfos:
            if pyfmo.orbitals.mulliken_contribution(sfo, mo) * 100 >= 1:
                sfo_energy = sfo.energy
                mo_energy = mo.energy
                energies.append([mo_energy, sfo_energy])

    return energies

def xvalues_contributions(energycontribution_y_f1, contribution_xy, contribution_xy_f1):

    result_list = []

    for pair in energycontribution_y_f1:
        found = False
        for x1 in contribution_xy:
            if pair[0] == x1[0]:
                for x2 in contribution_xy_f1:
                    if pair[1] == x2[0]:
                        result_list.append([x1[1]+0.08, x2[1]-0.08])
                        found = True
                        break
                if found:
                    break

    return result_list

def percentages(mos, sfos):

    percentages = []
    for mo in mos:
        for sfo in sfos:
            percentage = pyfmo.orbitals.mulliken_contribution(sfo, mo) * 100
            if percentage >= 1:
                percentages.append(percentage)

    return percentages

def population(sfos):

    populations = []
    for sfo in sfos:
        population = sum([pyfmo.Orbitals.mulliken_contribution(mo, sfo) * mo.occupation for mo in orbitals.mos])
        populations.append(population)

    return populations

def diagram(mos, orbitals, sfos_f1, sfos_f2, color='b'):

    # Creating y-values
    # MO
    energies_mo = energy_yvalues(mos)
    pair_occupations_mo = pair_occupation(mos)
    single_occupations_mo = single_occupation(mos)

    # SFO1
    energies_sfos_f1 = energy_yvalues(sfos_f1)
    pair_occupations_sfos_f1 = pair_occupation(sfos_f1)
    single_occupations_sfos_f1 = single_occupation(sfos_f1)

    # SFO2
    energies_sfos_f2 = energy_yvalues(sfos_f2)
    pair_occupations_sfos_f2 = pair_occupation(sfos_f2)
    single_occupations_sfos_f2 = single_occupation(sfos_f2)


    # Creating x-values
    # MO
    x1 = degenerate_xvalues(energies_mo)
    x2 = [(x + 1) for x in degenerate_xvalues(energies_sfos_f1)] 
    x3 = [(x - 1) for x in degenerate_xvalues(energies_sfos_f2)]

    # Pair SFO1 and SFO2
    x4 = occupation_xvalues(pair_occupations_mo)
    x5 = [(x + 1) for x in occupation_xvalues(pair_occupations_sfos_f1)]
    x6 = [(x - 1) for x in occupation_xvalues(pair_occupations_sfos_f2)]

    # Single SFO1 and SFO2
    x7 = occupation_xvalues(single_occupations_mo)
    x8 = [(x + 1) for x in occupation_xvalues(single_occupations_sfos_f1)]
    x9 = [(x - 1) for x in occupation_xvalues(single_occupations_sfos_f2)]


    # Plotting the scatters in the same plot
    ax1 = plt.gca()

    # MO
    ax1.scatter(x1, energies_mo, s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x4, pair_occupations_mo, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if single_occupations_mo != []:
        x1.scatter(x7, single_occupations_mo, s=200, marker="$↑$", linewidth=0.2, c='k')

    # SFO1
    ax1.scatter(x2, energies_sfos_f1, s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x5, pair_occupations_sfos_f1, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if single_occupations_sfos_f1 != []:
        ax1.scatter(x8, single_occupations_sfos_f1, s=200, marker="$↑$", linewidth=0.2, c='k')

    # SFO2
    ax1.scatter(x3, energies_sfos_f2,  s=1444, marker="_", linewidth=3, c='k')
    ax1.scatter(x6, pair_occupations_sfos_f2, s=200, marker="$⇅$", linewidth=0.2, c='k')
    if single_occupations_sfos_f2 != []:
        ax1.scatter(x9, single_occupations_sfos_f2, s=200, marker="$↑$", linewidth=0.5, c='k')

    # Plotting the lines
    contribution_xy = [list(x) for x in zip(energies_mo, x1)]
    contribution_xy_f1 = [list(x) for x in zip(energies_sfos_f1, x2)]
    contribution_xy_f2 = [list(x) for x in zip(energies_sfos_f2, x3)]

    # generate the x and y-values
    energycontribution_y_f1 = energies_contribution(mos, sfos_f1)
    energycontribution_y_f2 = energies_contribution(sfos_f2, mos)
    energycontribution_x_f1 = xvalues_contributions(energycontribution_y_f1, contribution_xy, contribution_xy_f1)
    energycontribution_x_f2 = xvalues_contributions(energycontribution_y_f2, contribution_xy_f2, contribution_xy)

    # Generate percentages to label the lines
    percentage_f1 = percentages(mos, sfos_f1)
    percentage_f2 = percentages(mos, sfos_f2)

    # Plotting the lines
    xvals = []
    for xpair, ypair, percentage in zip(energycontribution_x_f1, energycontribution_y_f1, percentage_f1):
        ax1.plot(xpair, ypair, alpha=0.5, ls='dashed', c=color, label=f"{int(percentage)/10 * 10} %")
        xvals.append(0.5)

    for xpair2, ypair2, percentage2 in zip(energycontribution_x_f2, energycontribution_y_f2, percentage_f2):
        ax1.plot(xpair2, ypair2, alpha=0.5, ls='dashed', c=color, label=f"{int(percentage2)/10 * 10} %")
        xvals.append(-0.5)

    # inline label for the contribution
    labelLines(ax1.get_lines(), align=True, xvals=xvals, ha="center", va="center", backgroundcolor="none", fontsize=6)


    # Getting the populations
    pop_f1 = population(sfos_f1)
    pop_f2 = population(sfos_f2)


    # Plotting the orbital label and energies for MO, SFO1 and SFO2
    # MO
    for i, j, l in zip(x1, energies_mo, mos):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(40,4), size=7, ha="center", va="top", textcoords="offset points")

    # SFO1
    for i, j, l, pop in zip(x2, energies_sfos_f1, sfos_f1, pop_f1):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(40,8), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")
        ax1.annotate(f"{round(pop, 2)} e", xy=(i, j), xytext=(40,-2), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")

    # SFO2
    for i, j, l, pop in zip(x3, energies_sfos_f2, sfos_f2, pop_f2):
        ax1.annotate(l, xy=(i, j), xytext=(0,-10), size=7, ha="center", va="top", textcoords="offset points")
        ax1.annotate(f"{round(j, 2)} eV", xy=(i, j), xytext=(-40,8), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")
        ax1.annotate(f"{round(pop, 2)} e", xy=(i, j), xytext=(-40,-2), size=7, ha="center", va="top", textcoords="offset points", backgroundcolor="w")

    #Tuning Axis
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.xticks([-1, 0, 1], ['SFO1', "MO", 'SFO2'])
    plt.ylabel(r'$\epsilon$ / eV')
    plt.xlim(-1.7, 1.7)



