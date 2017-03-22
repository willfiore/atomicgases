#include <iostream>  // input / output from console
#include <vector>    // variable sized containers to store flip times and states
#include <algorithm> // std::min, max
#include <cmath>     // math library
#include <numeric>   // for std::accumulate (summing over containers)
#include <fstream>   // file stream for saving to files
#include <sstream>   // string stream for creating the file name

#include <thread> // Multithreading
#include <mutex> // Locking vector while pushing

#include "Random.hpp"
#include "State.hpp"
#include "Plot.hpp"

std::mutex state_push_mutex;

// void generate_rates(state, rates)
// Given the current state vector of the system, update the corresponding rates
void generate_rates(std::vector<bool>& state, std::vector<double>& rates)
{
    const double pow_r_12 = pow(State::R, 12);

    // For each atom
    for (int k = 0; k < State::num_atoms; ++k) {

        // start interaction sum on 0
        double interaction_sum = 0;

        // For each atom
        for (int j = 0; j < State::num_atoms; ++j) {
            // Ignore the same atom (continue goes to beginning of loop again)
            if (j == k) continue;

            // Distance between two atoms is the shortest between real distance and wrapped distance
            // (to create periodic boundaries)
            int dist = min(std::abs(j - k), State::num_atoms - std::abs(j-k));

            // add to the interaction sum as given in the notes
            interaction_sum += double(state[j]) / double(pow(dist, 6));
        }

        // set the rate of this atom as given in the notes
        rates[k] = 1.f / (1.f + (pow_r_12)*pow(interaction_sum, 2));
    }
}

// double get_jump_time(rates)
// Get the next atomic jump time (from 0)
double get_jump_time(std::vector<double>& rates)
{
    // note: the std::accumulate method sums over the iterators given, starting from the value given
    // as a parameter ( std::accumulate(starting_iterator, ending_iterator, starting_value)
    return -log(Random::randomDouble(0, 1)) / std::accumulate(rates.begin(), rates.end(), 0.f);
}

// int get_jump_atom(rates)
// Get which atom performs the next jump
int get_jump_atom(std::vector<double>& rates)
{
    // First, sum over all the rates
    double sum_rates = std::accumulate(rates.begin(), rates.end(), 0.f);

    // Create a fixed-sized array of cumulatively summed rates
    std::vector<double> cum_rates(State::num_atoms, 0);

    // For each atom
    for (size_t i = 0; i < State::num_atoms; ++i) {
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
		std::vector<bool> current_state(State::num_atoms, false);

		// bool can either be true or false.
		// Here, true represents an atom in the Rydberg state,
		// false represents the ground state.
		// Converting bool to a numeric value (eg int, double) returns 0 or 1

		// Create a new array to store the transition rates of each atom, again
		// num_atoms in length
		std::vector<double> rates(State::num_atoms, 0);

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

		while (current_time < State::duration) {

			current_time += get_jump_time(rates);
			times.push_back(current_time);

			int flipped_atom = get_jump_atom(rates);
			// ! means not. Here we are setting the current state of the flipped atom
			// to the opposite of what it was previously.
			current_state[flipped_atom] = !current_state[flipped_atom];
			// push the state to the current_state vector.
			states.push_back(current_state);

			// Regenerate the rates
			generate_rates(current_state, rates);

			// ... repeat
		}

		// We have multiple threads trying to access global vectors in the State:: namespace.
		// This can cause problems if two threads try to access the same vector at the same time.
		// The lock guard makes the thread running this function own the mutex defined at the top
		// of the file. While the mutex is owned, no thread is able to proceed past this point
		// in its execution. The lock guard is released at the end of the block, after the
		// vectors have been pushed to.
		std::lock_guard<std::mutex> guard(state_push_mutex);

		std::cout << "C" << core << ": Simulated " << current_time << " seconds. (repeat " << State::current_repeat++ << ")" << std::endl;
		State::repeated_times.push_back(times);
		State::repeated_states.push_back(states);
	}
}

// Main entry-point for a C++ program
// This is where execution begins
int main()
{
    // Get some user input
    std::cout << "Number of atoms \t> ";
    std::cin >> State::num_atoms;

    std::cout << "Interaction range, R \t> ";
    std::cin >> State::R;

    std::cout << "Sim Duration (seconds) \t> ";
    std::cin >> State::duration;

    std::cout << "Num repeats \t> ";
    std::cin >> State::num_repeats;

	// Round threads to the next multiple of 4 (for nice division between threads)
	State::num_repeats = ((3 + State::num_repeats) / 4) * 4;

	// Split execution across 4 threads. Assuming a 4 core machine, this should
	// quadruple execution time.
	std::thread t1(generateData, 0, State::num_repeats/4);
	std::thread t2(generateData, 1, State::num_repeats/4);
	std::thread t3(generateData, 2, State::num_repeats/4);
	std::thread t4(generateData, 3, State::num_repeats/4);

	// Once all the data has generated, join the threads back together.
	// std::thread::join() pauses execution until the thread has finished.
	t1.join();
	t2.join();
	t3.join();
	t4.join();

	// Here, all threads have finished generating data, so we are safe to plot.

    // Plot
    Plot::init();

    Plot::plotStateGraph();
	Plot::newPlotWindow();
	Plot::plotDensityGraph();
	Plot::newPlotWindow();
	Plot::plotFluctuationGraph();
	Plot::newPlotWindow();
	Plot::plotSpatialCorrelations();
	Plot::newPlotWindow();
    Plot::plotAllSpatialCorrelations();
	Plot::newPlotWindow();

    // Pause so it doesn't exit immediately and we have time to see the graphs.
    system("pause");

    // The main function must return an integer to tell the operating system
    // whether the program ran successfully. In C++, a return value of 0
    // means the program exited successfully.
    return 0;
}