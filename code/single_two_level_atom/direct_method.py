# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:59:51 2017

@author: wfior
"""

import math
import random
import matplotlib.pyplot as plot
import numpy as np

k = 1

def time_from_probability (prob) :
    return ((-1/k) * math.log(1 - prob))
    
def generate_one_atom(n, p):
    times = [0]
    flips = [0]
    for i in range(n):
        time = time_from_probability(random.random()/2)
        times.append(time)
        times.append(time)
        
        flips.append(flips[-1])
        flips.append(1 if flips[-1] == 0 else 0)
    
    times.sort()
    
    if (p):
        plot.plot(times, flips)
        plot.ylim(-0.5, 1.5)
        
    return times, flips

def generate_curve(n, bins):
    times = []
    for i in range(n):
        times.extend(generate_one_atom(100, False)[0])
        
    min_time = 0
    max_time = max(times)
    
    final_times = np.linspace(min_time, max_time, bins)
    final_counts = []

    for i in range(bins-1):
        bin_min = final_times[i]
        bin_max = final_times[i+1]
        count = sum(bin_min < t < bin_max for t in times)
        final_counts.append(count)
        
    final_counts.append(0)
    
    final_counts = [x / n for x in final_counts]
    print(len(final_times))
    print(len(final_counts))
    
    plot.plot(final_times, final_counts)
        
generate_curve(20, 100)