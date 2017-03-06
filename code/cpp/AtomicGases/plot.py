# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 17:19:35 2017

@author: wfior
"""


import math
import csv
import matplotlib.pyplot as plt
import copy

data_og = []

with open('data.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        data_og.append(row)
        
data = copy.copy(data_og)

(num_repeats, R, num_atoms, duration) = map(int, data.pop(0))

repeated_densities = []

times = []
atoms = []



for repeat in range(num_repeats):
    times.append(data.pop(0))
    for atom in range(num_atoms):
        atoms.append(data.pop(0))
        
    # Calculate densities
    final_density = []
    for t in range(len(times)):
        density = 0
        for a in atoms:
            density += int(a[t])
            density /= len(atoms)
            
            final_density.append(density)
            
    repeated_densities.append(final_density)
        
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