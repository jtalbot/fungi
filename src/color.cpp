
#include "color.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

Spectrum Spectrum::loadSPD(std::string const& filename) {
    Spectrum result(0);

    std::ifstream f;
    f.open(filename, std::ifstream::in);

    unsigned int i = Spectrum::shortestWavelength;

    std::string line;
    for (std::string line; std::getline(f, line);) {
        float wavelength, weight;

        std::istringstream in(line);
        in >> wavelength >> weight;

        for (; i < Spectrum::longestWavelength && float(i) < wavelength; ++i) {
            result.spectrum[i - Spectrum::shortestWavelength] = weight;
        }
    }

    return result;
}
