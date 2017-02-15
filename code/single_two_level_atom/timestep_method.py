# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math
import random
import matplotlib.pyplot as plt

k = 1

def get_probability(time):
    return 0.5*(1-math.exp(-2*time*k))

def generate_state(plot=False):
    
    t = 0
    max_t = 3
    dt = 0.01

    time = [0]    
    state = [0]
    
    while t < max_t:
        t += dt
        
        if (plot):
            time.append(t)
            state.append(state[-1])
            
        time.append(t)
        state.append(1 if random.random() < get_probability(t) else 0)
        
    if (plot):
        plt.figure(figsize=(15, 5), dpi=80)
        plt.plot(time, state)
    return time, state

def average_state():
    N = 10000
    time = generate_state()[0]
    state = generate_state()[1]
    
    for i in range(N - 1):
        new_state = generate_state()[1]
        state = [sum(x) for x in zip(state, new_state)]
        
    state = [x / N for x in state]
            
    plt.figure(2)
    plt.plot(time, state)
    
    average_state = [get_probability(t) for t in time]
    plt.plot(time, average_state, 'r')
    
generate_state(True)
plt.ylim(-0.5, 1.5)
average_state()