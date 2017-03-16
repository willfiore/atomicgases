# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 17:19:35 2017

@author: wfior
"""


import math
import csv
import matplotlib.pyplot as plt
import numpy as np
import copy
import bisect

data_og = []

with open('data_100_8_1000000_1.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        data_og.append(row)
        
data = copy.copy(data_og)

(num_repeats, R, num_atoms, duration) = map(int, data.pop(0))

times = list(np.linspace(0.0001, duration, 1000))
atoms = []

for r in range(num_repeats):
    jump_times = []
    interp_atoms = []

    jump_times.extend(map(float, data.pop(0)))

    for a in range(num_atoms):
        extracted_atom = list(map(int, data.pop(0)))

        a_interp = [extracted_atom[bisect.bisect_left(jump_times, t)] for t in times]
        
        interp_atoms.append(a_interp)
    
    atoms.append(interp_atoms)
    
t_interp = np.linspace(0, duration, 1000)
a_interp

def plot_density():
    
    densities = [sum(x)/num_atoms for x in list(zip(*atoms[r])) for r in range(num_repeats)]
                 
    log_times = [math.log(t) for t in times]
    log_densities = [math.log(d) for d in densities]
    plt.plot(log_times, log_densities)
    
    # Takes first repeat only
def plot_excitation_graph():
    times = []
    atoms = []

    times.extend(map(float, data.pop(0)))
    for a in range(num_atoms):
        atoms.append(list(map(int, data.pop(0))))
        
        times_to_plot = [times[index] for (index, value) in enumerate(atoms[a]) if value == 1]
        atom_number = [a for x in times_to_plot]

        plt.scatter(times_to_plot, atom_number, color='black', marker=',', lw=0, s=1)
        
def plot_spatial_correlations():
    distances = list(range(1, num_atoms))
    
    g = []
    for x in distances:
        correlations = []
        n1s = []
        n1_xs = []

        for r in range(num_repeats):
            time = 50
            n1 = atoms[r][0][time]
            n1_x = atoms[r][x][time]

            print(n1_x)

            correlations.append(n1 * n1_x)
            n1s.append(n1)
            n1_xs.append(n1_x)
        
        average_correlation = sum(correlations) / num_repeats
        average_n1 = sum(n1s) / num_repeats
        average_n1_x = sum(n1_xs) / num_repeats

        g.append(average_correlation - (average_n1 * average_n1_x) )
        
    print(g)
    plt.plot(distances, g)
#         
#    plt.plot(x, g)
        
plot_density()
#plot_excitation_graph()
#plot_spatial_correlations()