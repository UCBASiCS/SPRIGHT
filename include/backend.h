#ifndef BACKEND_H
#define BACKEND_H

#include "binprocessor.h"
#include "frontend.h"
#include "step.h"
#include "utils.h"

#include <fftw3.h>
#include <set>
#include <unordered_map>

class BackEnd: public Step
{
private:
    const FrontEnd* frontEnd;
    BinProcessor* binProcessor;
    ffast_complex** observationMatrix;

    std::unordered_map<int,ffast_complex> decodedFrequencies;
    std::vector<double> realFrequenciesIndices;

    bool* changed;
    int maxBinUsed;
    std::set<int> realUsedSamples;

public:
    BackEnd(Chrono* newChrono, const Config* newConfig, const FrontEnd* newFrontEnd);
    ~BackEnd();
    void process();
    const std::unordered_map<int,ffast_complex>& getDecodedFrequencies() const;
    void swapDecodedFrequencies(std::unordered_map<int,ffast_complex>);
    const std::vector<double>& getRealFrequenciesIndices() const;
    void printObservations();
    int getMaxBinUsed();
    int getRealUsedSamplesNb();

private:
    void initialize();
    void peelFrom(int spectralIndex, ffast_complex amplitude, int binIndex);
    void getClusteredFrequencies();
    int getHash(int stage, int spectralIndex);
};

#endif // BACKEND_H
