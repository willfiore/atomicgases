#include <iostream>  // input / output from console
#include <array>     // fixed sized containers to store current states and rates
#include <vector>    // variable sized containers to store flip times and states
#include <cmath>     // math library
#include <numeric>   // for std::accumulate (summing over containers)
#include <random>    // random number generation
#include <fstream>   // file stream for saving to files

// Global variables 
const int num_atoms = 500;
const float duration = 100;

// These are provided by user input when the application is run
int R = 20;
int num_repeats = 1;

// float randomFloat(float a, float b)
// Returns a random floating-point in the interval [a, b]
float randomFloat(float a, float b)
{
    // create and seed a new random generator using the OS's random device
    std::mt19937 gen = std::mt19937(std::random_device()());
    // create a uniform real number distribution in the interval [a, b]
    std::uniform_real_distribution<float> dis(a, b);
    // use the random generator to produce a number from the distribution, and return it
    return dis(gen);
}

// void generate_rates(state, rates)
// Given the current state vector of the system, update the corresponding rates
void generate_rates(std::array<bool, num_atoms>& state, std::array<float, num_atoms>& rates)
{
    // For each atom
    for (int k = 0; k < num_atoms; ++k) {

        // start interaction sum on 0
        float interaction_sum = 0;

        // For each atom
        for (int j = 0; j < num_atoms; ++j) {
            // Ignore the same atom (continue goes to beginning of loop again)
            if (j == k) continue;

            // Distance between two atoms is the shortest between real distance and wrapped distance
            // (to create periodic boundaries)
            int dist = std::min(std::abs(j - k), num_atoms - std::abs(j-k));

            // add to the interaction sum as given in the notes
            interaction_sum += float(state[j]) / float(pow(dist, 6));
        }

        // set the rate of this atom as given in the notes
        rates[k] = 1.f / (1.f + (pow(R, 12))*pow(interaction_sum, 2));
    }
}

// float get_jump_time(rates)
// Get the next atomic jump time (from 0)
float get_jump_time(std::array<float, num_atoms>& rates)
{
    // note: the std::accumulate method sums over the iterators given, starting from the value given
    // as a parameter ( std::accumulate(starting_iterator, ending_iterator, starting_value)
    return -log(randomFloat(0, 1)) / std::accumulate(rates.begin(), rates.end(), 0.f);
}

// int get_jump_atom(rates)
// Get which atom performs the next jump
int get_jump_atom(std::array<float, num_atoms>& rates)
{
    // First, sum over all the rates
    float sum_rates = std::accumulate(rates.begin(), rates.end(), 0.f);

    // Create a fixed-sized array of cumulatively summed rates
    std::array<float, num_atoms> cum_rates;

    // For each atom
    for (size_t i = 0; i < num_atoms; ++i) {
        // Start the sum on 0
        float sum = 0;

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
    float r = randomFloat(0, 1);

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

// Main entry-point for a C++ program
// This is where execution begins
int main()
{
    // Get some user input
    std::cout << "R > ";
    std::cin >> R;
    std::cout << std::endl;

    std::cout << "Num repeats > ";
    std::cin >> num_repeats;
    std::cout << std::endl;

    // We want to output our data as a comma-seperated values (CSV) file,
    // to be read into python.
    std::ofstream file;    // Create a new file stream
    file.open("data.csv"); // Open the file

    // The operator << is used to push data into the file.

    // On the first line, push important information about the data-set
    // Note: the special character '\n' is a carriage return (new line)
    file << num_repeats << "," << R << "," << num_atoms << "," << duration << "\n";

    // Repeat num_repeats amount of times
    for (size_t r = 0; r < num_repeats; ++r) {
        // Print repeat number to console
        std::cout << "Repeat number " << r << std::endl;

        // Create a new array to store the current state.
        // std::array takes two template parameters as follows:
        // std::array<value_type, number_of_values>

        // 'bool' represents a boolean value, either true or false.
        // if 'bool' is cast to a numeric value like 'int' or 'float', then
        // true -> 1, false -> 0
        std::array<bool, num_atoms> current_state;

        // Here, true represents an atom in the Rydberg state,
        // false represents the ground state.

        // Initialize all atoms in the ground state to begin
        // (maybe there's a better way of doing this, I'm not sure. Who cares)
        // The syntax here means 'for every state in current_state'
        for (auto& state : current_state) {
            state = false;
        }

        // Create a new array to store the transition rates of each atom, again
        // num_atoms in length
        std::array<float, num_atoms> rates;

        // Generate the rates according to the states
        generate_rates(current_state, rates);

        // Create an empty vector of floats to store the jump times.
        std::vector<float> times;
        // We want to store the entire state of the system at each jump time, so
        // create another empty vector to store state arrays
        std::vector<std::array<bool, num_atoms> > states;

        // Start the time on 0
        // 'float' is a floating-point number, i.e one with decimal places
        // (as opposed to an integer)
        float current_time = 0;

        while (current_time < duration) {
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
            
            // Print percentage done
            std::cout << (100 * (current_time / duration)) << "%" << std::endl;

            // ... repeat
        }

        // Here, we are finished generating the data
        // Now, deal with saving it to the file:

        // For every time in times, push that value to the file, followed by a comma
        // We don't want the last value to have a comma after it so we only do this for
        // all times except the last one
        for (size_t t = 0; t < times.size() - 1; ++t) {
            file << times[t] << ",";
        }

        // Push the last time to the file, without the comma this time
        file << times[times.size() - 1];
        file << "\n"; // new line

        // For every atom
        for (size_t a = 0; a < num_atoms; ++a) {
            // For every time
            for (size_t t = 0; t < times.size() - 1; ++t) {
                // Push the state of atom a at time t to the file
                file << states[t][a] << ",";
            }

            // Again, push the last state without a comma
            file << states[times.size() - 1][a];
            file << "\n"; // new line
        }

        // Recall, all of this is inside a loop for the number of repeats
        // So, go back and do it all again. In the end, our file will look like:

        // - preliminary data
        // repeat 1: jump times
        // repeat 1: atom 0 states at each time
        // repeat 1: atom 1 states at each time
        // ...
        // repeat 1: atom N states at each time
        // repeat 2: jump times
        // repeat 2: atom 0 states at each time
        // ...
        // ...
        // repeat N: jump times
        // repeat N: atom jump times...
    }
    
    // Finally finished everything. Close the file stream (this saves the file, too)
    file.close();

    // The main function must return an integer to tell the operating system
    // whether the program ran successfully. In C++, a return value of 0
    // means the program exited successfully.
    return 0;
}