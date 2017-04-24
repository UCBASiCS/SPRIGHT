#ifndef FRONTEND_H
#define FRONTEND_H

#include "input.h"
#include "step.h"
#include "utils.h"

#include <set>
#include <vector>

class FrontEnd: public Step
{
private:
    std::unordered_map<int,ffast_complex> signal;

    ffast_complex** observationMatrix;
    int signalLength;

    bool countSamplesDone;  // Flag to indicate if the number of samples used by FFAST are counted.
    std::vector<int> delays;

public:
    FrontEnd(Chrono* newChrono, const Config* newConfig, const Input* newInput);
    ~FrontEnd();
    void process();
    ffast_complex** getObservationMatrix() const;
    int getUsedSamplesNb() const;
    std::vector<int> getDelays() const;
    const Input* input;
    std::set<int> usedSamples; // A set to keep track of samples used by FFAST

    const std::unordered_map<int,ffast_complex>& getTimeSignal() const;

private:
    void computeDelays();
};

#endif // FRONTEND_H
