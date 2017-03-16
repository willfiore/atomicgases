#ifndef STATE_H
#define STATE_H

#include <vector>

typedef std::vector<double> Times;
typedef std::vector<std::vector<bool> > States;

namespace State
{
    extern int num_atoms;
    extern int R;
    extern double duration;
    extern double real_duration;
    extern int num_repeats;

    extern std::vector<Times> repeated_times;
    extern std::vector<States> repeated_states;
};

#endif