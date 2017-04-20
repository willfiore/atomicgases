#ifndef RANDOM_H
#define RANDOM_H

#include <random>    // random number generation

class Random
{
public:
    static void seed(unsigned int);
    static double randomDouble(double a, double b);
	static bool randomBool();

private:
    static std::mt19937 generator;
};

#endif