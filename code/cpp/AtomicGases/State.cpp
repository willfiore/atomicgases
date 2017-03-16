#include "State.hpp"

namespace State
{
    int num_atoms = 500;
    int R = 20;
    double duration = 0;
    double real_duration = 0;
    int num_repeats = 1;
    int graph_type = 0;

    std::vector<Times> repeated_times;
    std::vector<States> repeated_states;
}