#ifndef PLOT_H
#define PLOT_H

#include "ext/gnuplot-iostream.h"

class Plot
{
public:
    static void init();
    static void plotStateGraph();

    static void plotDensityGraph();

    static void plotSpatialCorrelations(); // do last time
    static void plotSpatialCorrelations(double time);
private:
    static Gnuplot gp;
};

#endif