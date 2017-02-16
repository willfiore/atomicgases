### TODO:
### this method is currently not working correctly.

########################
### direct_method.py ###
########################

### Description:
# simulating transition probability for a single two-level atom

import math
import random
import matplotlib.pyplot as plot
import numpy as np

k = 1 # transition rate

def time_from_probability (prob):
    """ Returns the time at which a generated probability would place the
        atom in an excited state.
        
    Args:
        prob (float): probability of which to find the time for.
            
    Returns:
        time (float)
            
    """
    
    return ((-1/k) * math.log(1 - prob))
    
def generate_one_atom(n, p):
    """ Generate times and directions of state flips for a single atom
    
    Args:
        n (int): number of state flips to generate
        p (bool): Whether to plot the final result
        
    Returns:
        (times[], flips[]), where times[] is an array of flip times and
        flips[] is the direction of the corresponding flips.
    
    """
    
    # Create initial arrays to store times and flips
    times = [0]
    flips = [0]

    for i in range(n):
        
        # Generate time from random probability roll.
        time = time_from_probability(random.random()/2)
        
        times.append(time)
        times.append(time)
        
        # Flip always opposes the direction of the previous flip
        flips.append(flips[-1])
        flips.append(1 if flips[-1] == 0 else 0)
        
        # Note how the times and flips are doubled up so the resulting plot
        # consists of vertical, direct state changes instaed of diagonal lines
        # (interpolating flips makes no sense)
    
    times.sort()
    
    if (p):
        plot.plot(times, flips)
        plot.ylim(-0.5, 1.5)
        
    return times, flips

def generate_curve(n, bins):
    """ Generate average state for many atoms
    
    Args:
        n (int): number of atoms to simulate
        bins (int): Number of bins to use for averaging flips from multiple atoms
    """
    
    times = []
    for i in range(n):
        times.extend(generate_one_atom(100, False)[0])
        
    min_time = 0
    max_time = max(times)
    
    final_times = np.linspace(min_time, max_time, bins)
    final_counts = []

    # sum the total number of flips in this bin
    for i in range(bins-1):
        bin_min = final_times[i]
        bin_max = final_times[i+1]
        count = sum(bin_min < t < bin_max for t in times)
        final_counts.append(count)

    # make sure final_counts is the same length as final_times        
    final_counts.append(0)    
    final_counts = [x / n for x in final_counts] # normalize
    
    plot.plot(final_times, final_counts)
        
# generate curve
generate_curve(200, 1000)