#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <numeric>
#include <time.h>
#include <random>
#include <fstream>
#include <algorithm>

const int num_atoms = 500;
int R = 20;
int num_repeats = 1;
const float duration = 100;

float randomFloat(float a, float b)
{
    std::mt19937 gen = std::mt19937(std::random_device()());
    std::uniform_real_distribution<float> dis(a, b);
    return dis(gen);
}

void generate_rates(std::array<bool, num_atoms>& state, std::array<float, num_atoms>& rates)
{
    for (int k = 0; k < num_atoms; ++k) {

        float interaction_sum = 0;

        for (int j = 0; j < num_atoms; ++j) {
            if (j == k) continue;

            // Periodic boundaries
            int dist = std::min(std::abs(j - k), num_atoms - std::abs(j-k));

            interaction_sum += float(state[j]) / float(pow(dist, 6));
        }

        rates[k] = 1.f / (1.f + (pow(R, 12))*pow(interaction_sum, 2));
    }
}

float get_jump_time(std::array<float, num_atoms>& rates)
{
    return -log(randomFloat(0, 1)) / std::accumulate(rates.begin(), rates.end(), 0.f);
}

int get_jump_atom(std::array<float, num_atoms>& rates)
{
    float sum_rates = std::accumulate(rates.begin(), rates.end(), 0.f);

    std::array<float, num_atoms> cum_rates;
    for (size_t i = 0; i < cum_rates.size(); ++i) {

        float sum = 0;
        for (size_t j = 0; j <= i; ++j) {
            sum += rates[j];
        }
        cum_rates[i] = sum / sum_rates;
    }

    float r = randomFloat(0, 1);

    int atom = 0;
    for (auto& rate : cum_rates) {
        if (rate > r) break;
        atom++;
    }

    return atom;
}

int main()
{
    std::cout << "R > ";
    std::cin >> R;
    std::cout << std::endl;

    std::cout << "Num repeats > ";
    std::cin >> num_repeats;
    std::cout << std::endl;

    // Output as CSV
    std::ofstream file;
    file.open("data.csv");

    file << num_repeats << "," << R << "," << num_atoms << "," << duration << "\n";

    for (size_t r = 0; r < num_repeats; ++r) {
        std::cout << "Repeat number " << r << std::endl;

        std::array<bool, num_atoms> current_state;
        for (auto& i : current_state) {
            i = false;
        }

        std::array<float, num_atoms> rates;

        generate_rates(current_state, rates);

        std::vector<float> times;
        std::vector<std::array<bool, num_atoms> > states;

        float current_time = 0;
        while (current_time < duration) {
            current_time += get_jump_time(rates);
            times.push_back(current_time);

            int flipped_atom = get_jump_atom(rates);
            current_state[flipped_atom] = !current_state[flipped_atom];
            states.push_back(current_state);

            generate_rates(current_state, rates);
            std::cout << (100 * (current_time / duration)) << "%" << std::endl;
        }

        // Times
        for (size_t t = 0; t < times.size() - 1; ++t) {
            file << times[t] << ",";
        }
        file << times[times.size() - 1];
        file << "\n";

        for (size_t a = 0; a < num_atoms; ++a) {
            for (size_t t = 0; t < times.size() - 1; ++t) {
                file << states[t][a] << ",";
            }
            file << states[times.size() - 1][a];
            file << "\n";
        }
    }

    file.close();

    return 0;
}