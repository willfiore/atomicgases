import math
import random
import bisect
import numpy as np
import copy as copy
import matplotlib.pyplot as plt

num_atoms = 500

def generate_rates(R, state, rate):
    for k in range(num_atoms):
        list_of_atoms = list(range(num_atoms))
        list_of_atoms.remove(k)
        
        interaction_sum = 0
        for j in list_of_atoms:
            interaction_sum += state[j]/(abs(j-k)**6)
            
        rate[k] = 1/(1 + (R**12) * (interaction_sum**2))
        
def get_jump_time(rate):
    return -math.log(random.random())/sum(rate)
    
def get_jump_atom(rate):
    cum_rate = np.cumsum(rate) / sum(rate)
    return bisect.bisect_left(cum_rate, random.random())
        
def get_density(R):
    state = [0 for x in range(num_atoms)]
    rate = [0 for x in range(num_atoms)]
    
    generate_rates(R, state, rate)
    
    duration = 100
    
    last_jump_time = 0
    times = []
    total_state = []
    
    while (last_jump_time < duration):
        last_jump_time += get_jump_time(rate)
        times.append(last_jump_time)
        
        flipped_atom = get_jump_atom(rate)
        state[flipped_atom] = 1 if state[flipped_atom] == 0 else 0
        total_state.append(copy.copy(state))
        
        generate_rates(R, state, rate)    
        print(last_jump_time)
        
    num_times = 150
    q_times  = [x*(duration/num_times) for x in range(num_times)]
    q_states = []
    for i in range(num_times):
        q_states.append(total_state[bisect.bisect_left(times, q_times[i])])
        
    # density
    return q_times, [sum(x)/num_atoms for x in q_states]

def plot_average_density_graph(R):
    times = []
    densities = []
    num_repeats = 1
    for i in range(num_repeats):
        ret_vals = get_density(R)
        times = ret_vals[0]
        densities.append(ret_vals[1])
        
    final_densities = []
    for t in range(len(times)):
        final_density = 0
        for r in range(num_repeats):
            final_density += densities[r][t]
            
        final_density /= num_repeats
        final_densities.append(final_density)
        
    #plot density
    times.pop(0)
    final_densities.pop(0)
    times = [math.log(x) for x in times]
    final_densities = [math.log(x) for x in final_densities]
    plt.plot(times, final_densities)

#for R in np.linspace(4, 20, 3):
plot_average_density_graph(10)
    
plt.figure()

# plot data
t = []
p = []

#for i in range(len(q_times)):
#    t.extend([q_times[i] for n in range(q_states[i].count(1))])
#    p.extend([index for (index, value) in enumerate(q_states[i]) if value == 1])
#    
#plt.scatter(t, p, color='black',marker=',',lw=0, s=1)
#plt.xlabel("Time, t / s")
#plt.ylabel("Atom")