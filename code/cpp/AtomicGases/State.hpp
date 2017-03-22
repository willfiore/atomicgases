#ifndef STATE_H
#define STATE_H

#include <vector>
#include <string>

typedef char StateType;
typedef std::vector<double> Times;
typedef std::vector<std::vector<StateType> > States;

namespace State
{
    extern int num_atoms;
    extern int R;
    extern double duration;
    extern double real_duration;
    extern int num_repeats;

	extern int current_repeat;

    extern std::vector<Times> repeated_times;
    extern std::vector<States> repeated_states;

	std::string getInfoString();
};

#endif