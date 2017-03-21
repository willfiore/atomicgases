#include "Plot.hpp"
#include "State.hpp"

#include <algorithm>
#include <iterator>
#include <numeric> // std::accumulate
#include <iostream>

// Initialize GNUPlot with location
Gnuplot Plot::gp("\"C:\\Program files\\gnuplot\\bin\\gnuplot.exe\"");

void Plot::init()
{
    gp << "set samples 10000\n";
	gp << "set term wxt 0\n";
}

void Plot::newPlotWindow()
{
	static int x = 0;
	gp << "set term wxt " << ++x << "\n";
}

void Plot::plotStateGraph()
{
    Times& times = State::repeated_times[0];
    States& states = State::repeated_states[0];

    std::vector<double> x; // Time
    std::vector<int> y;   // Atom number

    // We only want to plot points for excited atoms:
    for (size_t t = 0; t < times.size(); ++t) {
        unsigned int total_excited = 0;
        
        for (size_t a = 0; a < State::num_atoms; ++a) {
            if (states[t][a]) {
                total_excited++;
                y.push_back(a);
            }
        }
        // Grow duplicated_times by number of excited atoms and fill with the jump time at that point
        x.insert(x.end(), total_excited, times[t]);
    }

	std::stringstream title;
	title << "Excited atoms in a 1D lattice over time. " << State::getInfoString();

	gp << "set title \"" << title.str() << "\"\n";
	gp << "set xlabel \"time, t / s\"\n";
	gp << "set ylabel \"Atom number\"\n";

    gp << "plot '-' with dots\n";
    gp.send1d(boost::make_tuple(x, y));
    gp.flush();
}

void Plot::plotDensityGraph()
{
	// Create evenly-spaced log times.
	std::vector<double> log_times;
	std::vector<double> log_density;

	int num_spaces = 2000;
	for (int i = 0; i < num_spaces; ++i) {
		log_times.push_back(double(i) / double(num_spaces) * log(State::duration));
	}

	// For every time in log_times
	for (auto& log_time : log_times) {

		double total_density = 0;

		// Undo the log time to find actual time
		double time = exp(log_time);

		// For every repeat
		for (size_t r = 0; r < State::num_repeats; ++r) {

			Times& times = State::repeated_times[r];
			States& states = State::repeated_states[r];

			// We want to find how many excited atoms are at this time..
			unsigned int total_excited = 0;

			// For every atom
			for (size_t a = 0; a < State::num_atoms; ++a) {

				// Get atom state at time 'time' (bisect.bisect_left equivalent)
				auto t_index = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), time));

				bool state;
				if (t_index >= times.size()) {
					state = states[states.size() - 1][a];
				}
				else {
					state = states[t_index][a];
				}

				// If it's excited, add to total excited state
				total_excited += state;
			}

			// Density is total excited divided by number of atoms
			double density = total_excited / double(State::num_atoms);

			// We want to average the density over many repeats so add this to a vector
			total_density += density;
		}

		double average_density = total_density / State::num_repeats;

		// We have found the average density corresponding to the value in log_times
		// Append it to log_density
		log_density.push_back(log(average_density));
	}

	std::stringstream title;
	title << "Density of excited atoms over time. " << State::getInfoString();

	gp << "set title \"" << title.str() << "\"\n";
	gp << "set xlabel \"time, t / s\"\n";
	gp << "set ylabel \"Density\"\n";

	gp << "plot '-' with lines\n";
	gp.send1d(boost::make_tuple(log_times, log_density));
	gp.flush();
}

void Plot::plotFluctuationGraph()
{
    // Create evenly-spaced times.
	std::vector<double> interp_times;
	std::vector<double> fluctuations;

    int num_spaces = 2000;
    for (int i = 0; i < num_spaces; ++i) {
        interp_times.push_back((double(i) / double(num_spaces)) * State::duration);
    }

    // For every time in interp_times
    for (auto& time : interp_times) {

        double total_density = 0;
		double total_sq_density = 0;

        // For every repeat
        for (size_t r = 0; r < State::num_repeats; ++r) {

            Times& jump_times = State::repeated_times[r];
            States& jump_states = State::repeated_states[r];

            // We want to find how many excited atoms are at this time..
            unsigned int total_excited = 0;

            // For every atom
            for (size_t a = 0; a < State::num_atoms; ++a) {

                // Get atom state at time 'time' (bisect.bisect_left equivalent)
                auto t_index = std::distance(jump_times.begin(), std::lower_bound(jump_times.begin(), jump_times.end(), time));

				bool state;
				if (t_index >= jump_times.size()) {
					state = jump_states[jump_states.size() - 1][a];
				}
				else {
					state = jump_states[t_index][a];
				}

                // If it's excited, add to total excited state
                total_excited += state;
            }

            // Density is total excited divided by number of atoms
            double density = total_excited / double(State::num_atoms);
			double sq_density = pow(density, 2);

            // We want to average the density over many repeats so add this to a vector
            total_density += density;
			total_sq_density += sq_density;
        }

        double average_density = total_density / State::num_repeats;
		double average_sq_density = total_sq_density / State::num_repeats;

        // We have found the average density corresponding to the value in log_times
        // Append it to log_density
		double fluctuation = sqrt(average_sq_density - pow(average_density, 2));
		fluctuations.push_back(fluctuation);
    }

	std::stringstream title;
	title << "Density fluctuation over time. " << State::getInfoString();

	gp << "set title \"" << title.str() << "\"\n";
	gp << "set xlabel \"time, t / s\"\n";
	gp << "set ylabel \"Density Fluctuation\"\n";

    gp << "plot '-' with lines\n";
    gp.send1d(boost::make_tuple(interp_times, fluctuations));
    gp.flush();
}

void Plot::plotSpatialCorrelations()
{
    plotSpatialCorrelations(State::duration);
}

void Plot::plotSpatialCorrelations(double time)
{
    std::vector<int> x;
    std::vector<double> g;

    // Need to calculate three expectation values:
    // 1) <n_1 * n_(1+x)>

	for (size_t i = 1; i < State::num_atoms / 2; ++i) {
		x.push_back(i);

		std::vector<int> correlations;
		std::vector<int> n_1;
		std::vector<int> n_x;

		for (size_t r = 0; r < State::num_repeats; ++r) {
			Times&  times = State::repeated_times[r];
			States& states = State::repeated_states[r];

			for (size_t s = 0; s < State::num_atoms; ++s) {

				// Find the staste at time 'time'
				auto t_index = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), time)) - 1;

				bool state_1 = State::repeated_states[r][t_index][s];
				bool state_x = State::repeated_states[r][t_index][(s + i) % State::num_atoms];

				correlations.push_back(state_1 * state_x);
				n_1.push_back(state_1);
				n_x.push_back(state_x);
			}
		}

		double average_correlation = std::accumulate(correlations.begin(), correlations.end(), 0.f) / (State::num_repeats * State::num_atoms);
		double average_n_1 = std::accumulate(n_1.begin(), n_1.end(), 0.f) / (State::num_repeats * State::num_atoms);
		double average_n_x = std::accumulate(n_x.begin(), n_x.end(), 0.f) / (State::num_repeats * State::num_atoms);

		double g_val = average_correlation - (average_n_1 * average_n_x);
		g.push_back(g_val);

	}

	std::stringstream title;
	title << "Spatial correlations at time " << time << ". " << State::getInfoString();

	gp << "set title \"" << title.str() << "\"\n";
	gp << "set xlabel \"time, t / s\"\n";
	gp << "set ylabel \"Density\"\n";

	gp << "plot '-' with lines\n";
	gp.send1d(boost::make_tuple(x, g));
	gp.flush();
}

void Plot::plotAllSpatialCorrelations()
{
	std::vector<double> interp_times;
	std::vector<int> x_vals;
	std::vector<std::vector<double> > g_vals;

	// Populating interp_times
	int num_spaces = 200;
	for (int i = 0; i < num_spaces; ++i) {
		double time = (double(i) / double(num_spaces)) * State::duration;
		interp_times.push_back(time);
		std::cout << time << std::endl;
	}

	// Populating x
	for (size_t i = 1; i < State::num_atoms / 2; ++i) {
		x_vals.push_back(i);
	}

	for (auto& time : interp_times) {
		std::vector<double> g;

		// Need to calculate three expectation values:
		// 1) <n_1 * n_(1+x)>

		for (size_t i = 1; i < State::num_atoms / 2; ++i) {

			std::vector<int> correlations;
			std::vector<int> n_1;
			std::vector<int> n_x;

			for (size_t r = 0; r < State::num_repeats; ++r) {
				Times&  times = State::repeated_times[r];
				States& states = State::repeated_states[r];

				for (size_t s = 0; s < State::num_atoms; ++s) {

					// Find the state at time 'time'
					auto t_index = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), time));

					bool state_1 = State::repeated_states[r][t_index][s];
					bool state_x = State::repeated_states[r][t_index][(s + i) % State::num_atoms];

					correlations.push_back(state_1 * state_x);
					n_1.push_back(state_1);
					n_x.push_back(state_x);
				}
			}

			double average_correlation = std::accumulate(correlations.begin(), correlations.end(), 0.f) / (State::num_repeats * State::num_atoms);
			double average_n_1 = std::accumulate(n_1.begin(), n_1.end(), 0.f) / (State::num_repeats * State::num_atoms);
			double average_n_x = std::accumulate(n_x.begin(), n_x.end(), 0.f) / (State::num_repeats * State::num_atoms);

			double g_val = average_correlation - (average_n_1 * average_n_x);
			g.push_back(g_val);
		}

		g_vals.push_back(g);
	}

	std::vector<std::vector<double> > plot_x(x_vals.size(), interp_times);
	std::vector<std::vector<int> > plot_y(interp_times.size(), x_vals);

	std::vector<std::vector<int> > plot_x_trans(interp_times.size(), std::vector<int>());

	for (size_t i = 0; i < plot_x.size(); ++i) {
		for (size_t j = 0; j < plot_x[i].size(); ++j) {
			plot_x_trans[j].push_back(plot_x[i][j]);
		}
	}

	std::stringstream title;
	title << "Spatial correlation map for varying time t and atom distance x. " << State::getInfoString();

	gp << "set title \"" << title.str() << "\"\n";
	gp << "set xlabel \"time, t / s\"\n";
	gp << "set ylabel \"Density\"\n";

	gp << "set view map\n";
	gp << "splot '-' with image\n";
	gp.send2d(boost::make_tuple(plot_x_trans, plot_y, g_vals));
	gp.flush();
}