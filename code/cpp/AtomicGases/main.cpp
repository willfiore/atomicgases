#include <iostream>  // input / output from console
#include <vector>    // variable sized containers to store flip times and states
#include <algorithm> // std::min, max
#include <cmath>     // math library
#include <numeric>   // for std::accumulate (summing over containers)
#include <fstream>   // file stream for saving to files
#include <sstream>   // string stream for creating the file name

#include <thread> // Multithreading
#include <mutex> // Locking vector while pushing

#include "ext/gnuplot-iostream.h" // Plotting

#include "Random.hpp" // Generating random numbers (using Mersenne-twister)

// GLOBALS ---------------------------------------------
// -----------------------------------------------------

typedef char StateType;
// state values are either 0 or 1 but C++ compresses
// std::vector<bool> into a bool per bit, slowing access
// so instead we use char, which increases memory usage 8-fold
// but should give performance gains when accessing vectors holding
// this type

typedef std::vector<double> Times;
typedef std::vector<std::vector<StateType> > States;

int num_atoms = 500;
int R = 20;
double decay = 0;
double duration = 0;
double real_duration = 0;
int num_repeats = 1;
int current_repeat = 0;

std::vector<Times> repeated_times;
std::vector<States> repeated_states;

std::mutex state_push_mutex;

// -----------------------------------------------------
// -----------------------------------------------------

// void generate_rates(state, rates)
// Given the current state vector of the system, update the corresponding rates
void generate_rates(std::vector<StateType>& state, std::vector<double>& rates)
{
    // No need to compute pow() multiple times. Store it as a constant
    const double pow_r_12 = pow(R, 12);

    // For each atom
    for (int k = 0; k < num_atoms; ++k) {

        // start interaction sum on 0
        double interaction_sum = 0;

        // For each atom
        for (int j = 0; j < num_atoms; ++j) {
            // Ignore the same atom
            if (j == k) continue;

            // Distance between two atoms is the shortest between real distance and wrapped distance
            // (to create periodic boundaries)
            int dist = min(std::abs(j - k), num_atoms - std::abs(j-k));

            // add to the interaction sum as given in the notes
            interaction_sum += double(state[j]) / double(pow(dist, 6));
        }

        // set the rate of this atom as given in the notes
        rates[k] = 1.f / (1.f + (pow_r_12)*pow(interaction_sum, 2));

        // account for decay by adding a linear decay value to 
        // the transition if the atom is in the excited state
		if (state[k] >= 1) {
			rates[k] += decay;
		}
    }
}

// double get_jump_time(rates)
// Get the next atomic jump time (from 0)
double get_jump_time(std::vector<double>& rates)
{
    return -log(Random::randomDouble(0, 1)) / std::accumulate(rates.begin(), rates.end(), 0.f);
}

// int get_jump_atom(rates)
// Get which atom performs the next jump
int get_jump_atom(std::vector<double>& rates)
{
    // First, sum over all the rates
    double sum_rates = std::accumulate(rates.begin(), rates.end(), 0.f);

    // Create a fixed-sized array of cumulatively summed rates
    std::vector<double> cum_rates(num_atoms, 0);

    // For each atom
    for (size_t i = 0; i < num_atoms; ++i) {
        // Start the sum on 0
        double sum = 0;

        // For each atom up to and including this atom
        for (size_t j = 0; j <= i; ++j) {
            // add the rate of that atom to the total sum
            sum += rates[j];
        }
        // Set the cumulative rate at this point in the array, and normalize it
        cum_rates[i] = sum / sum_rates;
    }

    // The cum_rates array now contains a series of values from 0 to 1
    // Generate a random number from 0 to 1
    double r = Random::randomDouble(0, 1);

    // Emulating bisect.bisect_left from Python:
    // Starting atom is the left-most atom in the array
    int atom = 0;
    // For every rate in cum_rates
    for (auto& rate : cum_rates) {
        // If r is inserted into the array, would this rate lie to the right of it?
        // If so, then we have found the bisection point. Break out of the loop.
        if (rate > r) break;

        // Otherwise, try the next value in cum_rates by advancing atom.
        atom++;
    }

    return atom;
}

// void generateData(core, repeats)
// Generate data "repeats" amount of times.
// 'Core' is a number from 0-3 indicating which CPU core
// is doing the processing. Just so it can be printed.
void generateData(int core, size_t repeats) {
	// Repeat num_repeats amount of times
	for (size_t r = 0; r < repeats; ++r) {

		// Create a new vector to store the current state.
		// std::vector takes one template parameter (between <>) indicating
		// which type of value it stores.

		// Here I have used the fill-constructor which has the following syntax:
		// std::vector<type> name (number_of_elements, default_value)
		std::vector<StateType> current_state(num_atoms, false);

		for (size_t i = 0; i < current_state.size(); i += R) {
            current_state[i] = true;
		}

		// Create a new array to store the transition rates of each atom, again
		// num_atoms in length
		std::vector<double> rates(num_atoms, 0);

		generate_rates(current_state, rates);

		// Create an empty vector of doubles to store the jump times.
		Times times;
		// We want to store the entire state of the system at each jump time, so
		// create another empty vector to store state arrays
		States states;

		// Start the time on 0
		// 'double' is a double precision floating-point number, (i.e one with decimal places)
		// (as opposed to an integer)
		double current_time = 0;

		// Insert first times and states
		times.push_back(0);
		states.push_back(current_state);

		time_t start_real_time, last_print_time, current_real_time;
		time(&start_real_time);
		time(&last_print_time);
        time(&current_real_time);

		double total_flips = 0;
		double excited_flips = 0;

		while (current_time < duration) {

            // 1) get next jump time
			current_time += get_jump_time(rates);
			times.push_back(current_time);

            // 2) get which atom jumps
			int flipped_atom = get_jump_atom(rates);
			current_state[flipped_atom] = !current_state[flipped_atom];
			total_flips+= 1;
			if (!current_state[flipped_atom]) excited_flips+= 1;
			// push the state to the current_state vector.
			states.push_back(current_state);

			// 3) regenerate the transition rates
			generate_rates(current_state, rates);

			// ... repeat
		}

		// We have multiple threads trying to access global vectors.
		// This can cause problems if two threads try to access the same vector at the same time.
		// The lock guard makes the thread running this function own the mutex defined at the top
		// of the file. While the mutex is owned, no thread is able to proceed past this point
		// in its execution. The lock guard is released at the end of the block, after the
		// vectors have been pushed to.
		std::lock_guard<std::mutex> guard(state_push_mutex);

		current_repeat++;
		int eta = difftime(current_real_time, start_real_time) * (num_repeats - current_repeat);

        // Print ETA
		std::cout << "C" << core << " in " << difftime(current_real_time, start_real_time) << ". " << num_repeats - current_repeat << " left." << std::endl;

        // Push data to final states
		repeated_times.push_back(times);
		repeated_states.push_back(states);
	}
}

// Main entry-point for a C++ program
// This is where execution begins
int main()
{
    // Get some user input
    std::cout << "Number of atoms \t> ";
    std::cin >> num_atoms;

    std::cout << "Interaction range, R \t> ";
    std::cin >> R;

    std::cout << "Sim Duration (seconds) \t> ";
    std::cin >> duration;

    std::cout << "Num repeats \t> ";
    std::cin >> num_repeats;

    double minDecay, maxDecay, decayIncrement;

    std::cout << "Min decay \t> ";
    std::cin >> minDecay;

    std::cout << "Max decay \t> ";
    std::cin >> maxDecay;

    std::cout << "Decay increment\t> ";
    std::cin >> decayIncrement;

    std::vector<double> decays;
    std::vector<double> stationary_state_densities;
    std::vector<double> stationary_state_fluctuations;

    // Set up gnuplot interface
    Gnuplot gp("\"C:\\Program files\\gnuplot\\bin\\gnuplot.exe\"");

    gp << "set samples 10000\n";
    gp << "set term wxt 0\n";

    int plot_window = 0;

    for (double x = minDecay; x <= maxDecay; x += decayIncrement) {

        // Clear previous data
        repeated_times.clear();
        repeated_states.clear();

        // Set correct decay
        decay = x;
        decays.push_back(decay);

        // Round threads to the next multiple of 4 (for nice division between threads)
        num_repeats = ((3 + num_repeats) / 4) * 4;

        // Split execution across 4 threads. Assuming a 4 core machine, this should
        // quadruple execution time.
        std::thread t1(generateData, 0, num_repeats / 4);
        std::thread t2(generateData, 1, num_repeats / 4);
        std::thread t3(generateData, 2, num_repeats / 4);
        std::thread t4(generateData, 3, num_repeats / 4);

        // Once all the data has generated, join the threads back together.
        // std::thread::join() pauses execution until the thread has finished.
        t1.join();
        t2.join();
        t3.join();
        t4.join();

        double average_density = 0;
        double average_sq_density = 0;

        // calculate density of stationary state
        for (size_t r = 0; r < num_repeats; ++r) {

            int total_excited = 0;

            for (size_t a = 0; a < num_atoms; ++a) {
                total_excited += repeated_states[r].back()[a];
            }

            double density = double(total_excited) / double(num_atoms);
            average_density += density;
            average_sq_density += pow(density, 2);
        }

        average_density /= double(num_repeats);
        average_sq_density /= double(num_repeats);

        stationary_state_densities.push_back(average_density);

        // calculate fluctuation of stationary state
        double fluctuation = (average_sq_density - pow(average_density, 2)) / average_density;
        stationary_state_fluctuations.push_back(fluctuation);
    }

    // Plot stationary state densities vs decay
    gp << "set term wxt " << ++plot_window << "title 'Stat state. density vs. decay rate'\n";
    gp << "plot '-' with lines\n";
    gp.send1d(boost::make_tuple(decays, stationary_state_densities));
    gp.flush();

    // Plot stationary state fluctuations vs decay
    gp << "set term wxt " << ++plot_window << "title 'Stat. state fluc. vs. decay rate'\n";
    gp << "plot '-' with lines\n";
    gp.send1d(boost::make_tuple(decays, stationary_state_fluctuations));
    gp.flush();

    // Pause so it doesn't exit immediately and we have time to see the graphs.
    system("pause");

    return 0;
}