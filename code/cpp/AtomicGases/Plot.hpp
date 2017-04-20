#ifndef PLOT_H
#define PLOT_H

#include "ext/gnuplot-iostream.h"

class Plot
{
public:
    static Gnuplot gp;

	static void init();
	static void newPlotWindow();

    static void plotStateGraph(int repeat);

    static void plotDensityGraph();
	static void plotFluctuationGraph();

    static void plotSpatialCorrelations(); // do last time
    static void plotSpatialCorrelations(double time);
	static void plotAllSpatialCorrelations();
};

#endif