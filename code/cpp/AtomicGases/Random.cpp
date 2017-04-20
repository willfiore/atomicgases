#include "Random.hpp"

// Initialize random generator with seed determined by random device
std::mt19937 Random::generator = std::mt19937(std::random_device()());

// Random::seed
// Re-seed the random generator with a specific seed
void Random::seed(unsigned int seed)
{
    generator.seed(seed);
}

// Random::randomFloat(float a, float b)
// Generate a random floating-point number in the interval [a, b]
double Random::randomDouble(double a, double b)
{
    // create a uniform real number distribution in the interval [a, b]
    std::uniform_real_distribution<double> distribution(a, b);
    // use the random generator to produce a number from the distribution, and return it
    return distribution(generator);
}

bool Random::randomBool()
{
	std::uniform_int_distribution<int> distribution(0, 1);
	return distribution(generator);
}