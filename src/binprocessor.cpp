#include "binprocessor.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>

BinProcessor::BinProcessor(Chrono *newChrono, const Config* newConfig, ffast_complex** newObservationMatrix, const std::vector<int>& newDelays):
    Step(newChrono, newConfig), observationMatrix(newObservationMatrix), delays(newDelays)
{
    signalLength     = config->getSignalLength();
    delaysNb         = config->getDelaysNb();
    chainsNb         = config->getChainsNb();
    delaysPerBunchNb = config->getDelaysPerBunchNb();

    MLdetection         = config->needToUseMaximumLikelihoodDetection();

    signalVector    = (ffast_complex*) malloc(delaysNb * sizeof(ffast_complex));
    thresholds      = (ffast_real*) malloc(config->getBinsNb() * sizeof(ffast_real));
    directionVector = (ffast_complex*) malloc(delaysNb * sizeof(ffast_complex));
    computeThresholds();
}

BinProcessor::~BinProcessor()
{
    free(signalVector);
    free(thresholds);
    free(directionVector);
}

void BinProcessor::process()
{
    if (MLdetection)
    {
        // MLprocess tries out all possible locations and estimates bin signal
        MLprocess();
    }
    else
    {
        computeLocation();
        estimateBinSignal();
    }
}

void BinProcessor::adjustTo(int newBinAbsoluteIndex, int newBinRelativeIndex, int newStage)
{
    binAbsoluteIndex = newBinAbsoluteIndex;
    binRelativeIndex = newBinRelativeIndex;
    stage   = newStage;
    samplingPattern = config->getSamplingPattern(stage);
}

bool BinProcessor::isSingleton()
{
    // a singleton is not zero-ton
    if (isZeroTon())
    {
        return false;
    }
    process();
    
    // if the bin is a singleton, when we peel the frequency from the bin
    // the remaining power should be coming from noise only
    // also we need the singleton to have energy larger than a threshold (minimumEnergy)

    if (noise / delaysNb <= thresholds[stage] && config->getDelaysNb()*std::norm(amplitude) > minimumEnergy)
    {
        // The bin is a singleton
        return true;
    }

    // The bin is a 'multi-ton'
    return false;
}

int BinProcessor::getLocation() const
{
    
    return location;
}

ffast_complex BinProcessor::getAmplitude() const
{
    return amplitude;
}

ffast_complex BinProcessor::getSignal(int delayIndex) const
{
    return signalVector[delayIndex];
}

bool BinProcessor::isZeroTon() const
{
    ffast_real energy = 0;

    for (int delayIndex=0; delayIndex<delaysNb; delayIndex++)
    {
        energy += std::norm(observationMatrix[binAbsoluteIndex][delayIndex]);
    }
    if (energy <= thresholds[stage])
    {
        return true;
    }

    return false;
}

// bool BinProcessor::isSingleton()

void BinProcessor::MLprocess()
{
    ffast_real MLnoise = std::numeric_limits<ffast_real>::infinity();
    int MLlocation = 0;

    /*
        Iterate over all balls that could map to this bin, and find the best one.
    */

    int homSolution = mapToInt(binRelativeIndex, config->getBinNumIndices(stage), samplingPattern, config->getSignalBitLength());
    int locGenerator = pow(2,(config->getSignalBitLength())) -1 - samplingPattern;

    int nnz = config->getSignalBitLength() - config->getBinNumIndices(stage);
    for (int i =0; i <pow(2,nnz); i++)
    {
        location = homSolution ^ mapToInt(i,nnz,locGenerator,config->getSignalBitLength());

        estimateBinSignal();
        if (noise < MLnoise)
        {
            MLnoise = noise;
            MLlocation = location;
        }
    }
    location = MLlocation;
    noise = MLnoise;
    estimateBinSignal();
    return ;
}

void BinProcessor::estimateBinSignal()
{
    int delayIndex = 0;
    amplitude      = 0;

    for(auto delayIterator=delays.begin(); delayIterator != delays.end(); ++delayIterator)
    {
        directionVector[delayIndex] = (numOnes(*delayIterator & location) % 2) ? -1 : 1;
        amplitude += directionVector[delayIndex] * observationMatrix[binAbsoluteIndex][delayIndex];
        delayIndex++;
    }

    // average the found amplitudes
    amplitude /= delaysNb;

    noise = 0;

    for (int delayIndex=0; delayIndex<delaysNb; delayIndex++)
    {
        // compute the bin signal that would result from a tone with such amplitude and location
        signalVector[delayIndex] = amplitude*directionVector[delayIndex];
        // the residual between the actual bin measurement and the computed one is the noise
        noise += std::norm(observationMatrix[binAbsoluteIndex][delayIndex]-signalVector[delayIndex]);
    } 
}

void BinProcessor::computeLocation()
{
    int tempLoc = 0;

    for(int q = 0; q < config->getSignalBitLength(); q++) {
        int runningSum = 0;
        for (int i = 0; i <config->getChainsNb(); i++) {
            ffast_complex Dp = observationMatrix[binAbsoluteIndex][(config->getSignalBitLength()+1)*i];
            ffast_complex Dpq = observationMatrix[binAbsoluteIndex][(config->getSignalBitLength()+1)*i+q+1];
            bool sgn1 = std::real(Dp) > 0;
            bool sgn2 = std::real(Dpq) > 0;
            bool sgn = sgn1 ^ sgn2;
            runningSum+= 2*sgn - 1;            
        }
        if (runningSum >=0)
            tempLoc += 1<<(config->getSignalBitLength() - 1- q);
    }
    location = tempLoc;
}

void BinProcessor::computeThresholds()
{
    // for computing the histogram for zeroton-singleton-multiton test
    std::vector<double> energyBins;
    double noiseEstimation = 0;
    double tempEnergy = 0;
    int    energyHistogramBinsCounted = 0;

    if ( config->isNoisy() )
    {
        // go over the stages

        for (int stage = 0; stage < config->getBinsNb(); ++stage)
        {
            // go over the bins in the stage
            for (int i = 0; i < config->getBinNumItems(stage); ++i)
            {
                // this is the energy of a bin measurement
                tempEnergy = 0;
                // go over the delays (go over bin measurements)
                for (int j = 0; j < config->getDelaysNb(); j++)
                {
                    // std::norm calculates norm-squared
                    tempEnergy += std::norm(observationMatrix[config->getBinOffset(stage)+i][j]);
                }
                energyBins.push_back(tempEnergy);
            }
        }

        // make histogram below
        // sort the energy bins
        std::sort (energyBins.begin(), energyBins.end());

        bool noiseLevelCrossed = false;
        while (!noiseLevelCrossed)
        {
            noiseEstimation += energyBins[energyHistogramBinsCounted];
            energyHistogramBinsCounted++;

            // if the bin energy is 10 times the previous one, declare it has signal in it
            if ( energyBins[energyHistogramBinsCounted] / energyBins[energyHistogramBinsCounted-1] >= 10 )
            {
                noiseLevelCrossed = true;
            }
        }



        // CHOICE 1: this corresponds to the average energy of the zero-ton bins
        noiseEstimation /= energyHistogramBinsCounted;

        // CHOICE 2: this corresponds to the maximum energy of the zero-ton bins
        // noiseEstimation = energyBins[energyHistogramBinsCounted-1];

    }

    // Minimum energy for the signal to be accepted as non-zero
    // 10 percent of the minimum SNR

    if( !config->isNoisy() ) // noiseless ongrid
    {
        minimumEnergy = pow(10,-8);
    }
    else if ( config->isNoisy() ) // noisy ongrid
    {
        /*
            We want a threshold on the minimum signal energy to eliminate false
            detections. A simple one would be of the sort:
            minimumEnergy = 0.1 * noiseEstimation * pow(10,config->getSNRdB()/10);
            However, for large signal to noise ratio, this is not working since there
            is always a noise floor around 1e-18 due to numerical issues, hence when
            multiplied with high SNR it blows up. Hence, for small SNR we do that, but
            for large SNR we clip it.
            We chose it to clip to 1000 times the noise floor.
        */
        minimumEnergy = std::min( 0.1 * noiseEstimation * pow(10,config->getSNRdB()/10) , 1000 * noiseEstimation );

    }

    /*
        Base threshold for the remaining bin signal to be considered as
        zero-ton after the signal is peeled from it. This value, in theory, is zero
        for noiseless simulations. However, due to machine precision, it needs to be positive
        value. For noiseless simulations we chose it to be 1e-13. If this value is high,
        we start to get false detections. One can decrease the thrsholds or increase the
        number of delays to reduce the false detections.
    */
    ffast_real baseThreshold = pow(10, -13);

    ffast_real factor;

    // Thresholds are obtained using tailbounds of Gaussian random variable
    if (delaysNb < 10)
    {
        factor = 4;
    }
    else if (delaysNb < 20)
    {
        factor = 3;
    }
    else if (delaysNb < 50)
    {
        factor = 2;
    }
    else
    {
        factor = 1.5;
    }

    for (int stage=0; stage<config->getBinsNb(); stage++)
    {
        thresholds[stage] = baseThreshold;

        if (config->isNoisy() )
            thresholds[stage] = pow(10,-10) + (ffast_real) factor * noiseEstimation;
    }
}
