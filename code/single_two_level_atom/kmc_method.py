# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math
import random
import matplotlib.pyplot as plt
import bisect

k = 1

def get_jump_time():
    return -(1/k)*math.log(random.random())
    
def simulate_atom(duration):
    time = [0]
    state = [0]
    
    while time[-1] < duration:
        time.append(get_jump_time() + time[-1])
        state.append(1 if state[-1] == 0 else 0)
        
    if (time[-1] > duration):
        time.pop()
        state.pop()
    
    return time, state
    
def get_average(N, duration, num_bins):
    average_state = [0 for x in range(num_bins)]
    time_x = [duration*(x/num_bins) for x in range(num_bins)]
    
    for a in range(N):
        single_atom_time, single_atom_state = simulate_atom(duration)
        
        for b in range(1,num_bins):
            index = bisect.bisect_left(single_atom_time, time_x[b]) - 1
            average_state[b] += single_atom_state[index]
            
    average_state = [a/N for a in average_state]
    
    print(average_state)
    plt.plot(time_x, average_state)

get_average(100000, 5, 1000)