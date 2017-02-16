##########################
### timestep_method.py ###
##########################

### Description:
# simulating transition probability for a single two-level atom

import math
import random
import matplotlib.pyplot as plt

k = 1 # transition rate

def get_excited_probability(time):
    """ Gets the probability of a single atom being in an excited state at
        a specific time.
    
    Args:
        time (float): time at which to find the probability
        
    """
    return 0.5*(1-math.exp(-2*time*k))

def generate_state(max_t=3, dt=0.01, plot=False):    
    """ Generates the state over time for a single atom
    
    Args:
        max_t (float): duration of time to simulate (in seconds)
        dt    (float): timestep interval (in seconds)
        plot   (bool): whether to plot the result afterwards
        
    Returns:
        (time[], state[]), where time[] is an array of all generated time values
        from 0 to max_t, and state[] is the corresponding states (either 0 or 1)
        for each time.
        
    """
    
    # starting time
    t = 0
    
    time = [t]
    state = [0]
    
    while t < max_t:
        t += dt
        
        # when plotting, we want to double up points so straight vertical lines
        # are drawn between transitions as opposed to diagonals (interpolating
        # between the two states 0 and 1 doesn't make sense)
        if (plot):
            time.append(t)
            state.append(state[-1])
            
        time.append(t)
        # generate a random probability between 0 and 1 and see if its below
        # the threshold for being in the excited state. then put the corresponding
        # state into the state array.
        state.append(1 if random.random() < get_excited_probability(t) else 0)
        
    if (plot):
        plt.figure(figsize=(8, 3), dpi=80)
        plt.plot(time, state)
        
    return time, state

def average_state(N=10000):
    """ Averages the state over many atoms and plots the result
    
    Args:
        N (int): number of atoms to simulate
    
    """
    
    # Grab initial time and state
    time = generate_state()[0]
    state = generate_state()[1]

    for i in range(N - 1):
        new_state = generate_state()[1]
        # Perform elementwise addition on the existing states with new states
        state = [sum(x) for x in zip(state, new_state)]
        
    # Normalise the final state
    state = [x / N for x in state]

    plt.figure(2)
    plt.plot(time, state)
    
    # Compare the plot with the exact probability function
    average_state = [get_excited_probability(t) for t in time]
    plt.plot(time, average_state, 'r')

# Generate single example state, and plot it
max_t = 2
generate_state(max_t, 0.01, True)
plt.ylim(-0.5, 1.5)
plt.xlim(0, max_t)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off

# Generate average state, and plot it
average_state(10000)