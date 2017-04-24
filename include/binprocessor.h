#ifndef BINPROCESSOR_H
#define BINPROCESSOR_H

#include "step.h"
#include "utils.h"

#include <vector>

class BinProcessor: public Step
{
private:
    ffast_complex** observationMatrix;
    int signalLength;
    int delaysNb;           
    int chainsNb;
    int delaysPerBunchNb;

    ffast_real minimumEnergy;

    int location;
    ffast_complex amplitude;
    ffast_complex* signalVector;

    int binAbsoluteIndex;
    int binRelativeIndex;
    int stage;
    int samplingPattern;
    bool MLdetection;
    ffast_real noise;
    ffast_real* thresholds;
    ffast_complex* directionVector;
    std::vector<int> delays;

public:
    BinProcessor(Chrono* newChrono, const Config* newConfig, ffast_complex** newObservationMatrix, const std::vector<int>& newDelays);
    ~BinProcessor();
    void process();
    void adjustTo(int newBinAbsoluteIndex, int newBinRelativeIndex, int newStage);
    bool isSingleton();
    int getLocation() const;
    ffast_complex getAmplitude() const;
    ffast_complex getSignal(int delayIndex) const;

private:
    bool isZeroTon() const;
    void MLprocess();
    void estimateBinSignal();
    void computeLocation();
    void computeThresholds();
};

#endif // BINPROCESSOR_H
