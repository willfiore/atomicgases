#include "State.hpp"

#include <sstream>

namespace State
{
    int num_atoms = 500;
    int R = 20;
	double decay = 0;
    double duration = 0;
    double real_duration = 0;
    int num_repeats = 1;
    int graph_type = 0;

	int current_repeat = 0;

    std::vector<Times> repeated_times;
    std::vector<States> repeated_states;

	std::string getInfoString()
	{
		std::stringstream ss;
		ss << num_atoms << " atoms, R = " << R << ", " << num_repeats << " repeats.";
		return ss.str();
	}
}