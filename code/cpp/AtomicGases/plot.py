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

data_og = []

with open('data_500_20_100000_3.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        data_og.append(row)
        
data = copy.copy(data_og)

(num_repeats, R, num_atoms, duration) = map(int, data.pop(0))

def plot_density():
    log_t_interp = np.linspace(math.log(0.0001), math.log(duration), 50)
    log_d_interp = []
    
    for r in range(num_repeats):
        times = []
        atoms = []
        times.extend(map(float, data.pop(0)))
        for a in range(num_atoms):
            atoms.append(list(map(int, data.pop(0))))
            
        density = [sum(x)/num_atoms for x in list(zip(*atoms))]
    
        log_times = [math.log(t) for t in times]
        log_density = [math.log(d) for d in density]
    
        log_d_interp.append(list(np.interp(log_t_interp, log_times, log_density)))
        
    zipped_d_interp = list(zip(*log_d_interp))
    average_d = [sum(d) for d in zipped_d_interp]
    
    plt.plot(log_t_interp, average_d)
    
times = []
atoms = []
    # Takes first repeat only
def plot_excitation_graph():
    times.extend(map(float, data.pop(0)))
    for a in range(num_atoms):
        atoms.append(list(map(int, data.pop(0))))
        
        times_to_plot = [times[index] for (index, value) in enumerate(atoms[a]) if value == 1]
        atom_number = [a for x in times_to_plot]

        plt.scatter(times_to_plot, atom_number, color='black', marker=',', lw=0, s=1)
        
plot_density()
#plot_excitation_graph()
        
# Plot excitation graph
#for a in range(len(atoms)):
#    times_to_plot = [times[index] for (index, value) in enumerate(atoms[a]) if value == '1']
#    atom_number = [a for x in times_to_plot]
#    
#    plt.scatter(times_to_plot, atom_number, color='black',marker=',',lw=0, s=1)
    
#plt.figure()
    
# Plot density graph
#final_density = []
#for t in range(len(times)):
#    density = 0
#    for a in atoms:
#        density += int(a[t])
#
#    density /= len(atoms)
#    final_density.append(density)
#    
#times = [math.log(float(t)) for t in times]
#final_density = [math.log(d) for d in final_density]
#    
#plt.plot(times, final_density)