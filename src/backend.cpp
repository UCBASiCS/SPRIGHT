#include "backend.h"
#include "helper.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>


BackEnd::BackEnd(Chrono* newChrono, const Config* newConfig, const FrontEnd* newFrontEnd):
    Step(newChrono, newConfig), frontEnd(newFrontEnd)
{
    observationMatrix = frontEnd->getObservationMatrix();
    changed = (bool*) malloc(config->getBinsSum() * sizeof(bool));
    maxBinUsed = 0;
}

BackEnd::~BackEnd()
{
    free(changed);
}



void BackEnd::process()
{
    chrono->start("BackEnd");
    binProcessor = new BinProcessor(chrono, config, observationMatrix, frontEnd->getDelays());
    int binAbsoluteIndex;
    bool singletonFound = true;
    initialize();
    maxBinUsed =0;
    while(singletonFound)
    {
        singletonFound = false;
        for (int stage = 0; stage < config->getBinsNb(); stage++)
        {
            binAbsoluteIndex = config->getBinOffset(stage);

            for (int binRelativeIndex = 0; binRelativeIndex < config->getBinNumItems(stage); binRelativeIndex++)
            {
                // Copies the values of bin index and related variables in binProcessor object.
                binProcessor->adjustTo(binAbsoluteIndex, binRelativeIndex, stage);
                // only check for singletons in the updated bins
                if (
                    changed[binAbsoluteIndex]
                    && binProcessor->isSingleton()
                    && decodedFrequencies.count(binProcessor->getLocation()) == 0
                   )
                {
                    singletonFound = true;
                    decodedFrequencies[binProcessor->getLocation()] = binProcessor->getAmplitude();
                    peelFrom(binProcessor->getLocation(), binProcessor->getAmplitude(), binAbsoluteIndex);
                    changed[binAbsoluteIndex] = false;
                }

                // mark this bin as processed
                // remark: this state can be changed if a ball is peeled from this bin
                changed[binAbsoluteIndex] = false;

                if ((int) decodedFrequencies.size() == config->getSignalSparsityPeeling())
                {
                    goto stopPeeling;
                }

                binAbsoluteIndex++;
            }
        }
    }


    stopPeeling:

    if ( config->applyWindow() )
    {
        getClusteredFrequencies();
    }
    else
    {
        for (auto it = decodedFrequencies.cbegin(); it != decodedFrequencies.cend(); ++it)
        {
            realFrequenciesIndices.push_back(it->first);
        }
    }

    delete binProcessor;

    chrono->stop("BackEnd");
}

int BackEnd::getMaxBinUsed() {
    return maxBinUsed + 1;
}

int BackEnd::getRealUsedSamplesNb() {
    return realUsedSamples.size();
}

void BackEnd::printObservations() {
    for (int stagePrint = 0; stagePrint < config->getBinsNb(); stagePrint++) {
        std::cout<< "-------------------------------" << std::endl;
        for (int r = 0; r < config->getBinNumItems(stagePrint); r++) {
            for (int c = 0; c < config->getDelaysNb(); c++) {
                std::cout << observationMatrix[r + config->getBinOffset(stagePrint)][c] << ", ";
            }
            std::cout << std::endl;
        }
    }
}

const std::vector<double>& BackEnd::getRealFrequenciesIndices() const
{
    return realFrequenciesIndices;
}

const std::unordered_map<int,ffast_complex>& BackEnd::getDecodedFrequencies() const
{
    return decodedFrequencies;
}

void BackEnd::swapDecodedFrequencies(std::unordered_map<int,ffast_complex> newDecodedFrequencies)
{
    decodedFrequencies.swap(newDecodedFrequencies);
}

void BackEnd::initialize()
{
    decodedFrequencies.clear();
    realFrequenciesIndices.clear();

    // peeling engine checks for singletons only in the updated (changed) bins
    // initialize all to true so that all the bins are check on the first run through
    for (int binAbsoluteIndex = 0; binAbsoluteIndex < config->getBinsSum(); binAbsoluteIndex++)
    {
        changed[binAbsoluteIndex] = true;
    }
}

void BackEnd::peelFrom(int spectralIndex, ffast_complex amplitude, int binIndex)
{
    int hash;
    // peel the contribution of the ball from all stages
    
    std::vector<int> delays = frontEnd->getDelays();
    ffast_complex* toSubtract = (ffast_complex *) malloc(config->getDelaysNb() * sizeof(ffast_complex));
    
    for (int i =0; i < config->getDelaysNb(); i++) {
        int currDelay = delays[i];
        toSubtract[i] = amplitude * (numOnes(currDelay & spectralIndex) % 2 ? -1.0 : 1.0);
    }

    for (int stage = 0; stage < config->getBinsNb(); stage++)
    {
        // index of the bin the ball is connected to at this stage
        hash = getHash(stage, spectralIndex);
        for (int delayIndex = 0; delayIndex < config->getDelaysNb(); delayIndex++)
        {
            observationMatrix[hash + config->getBinOffset(stage)][delayIndex] -= toSubtract[delayIndex];
            if (hash + config->getBinOffset(stage) == binIndex)
                observationMatrix[hash + config->getBinOffset(stage)][delayIndex] = 0;
            
        }
        changed[hash + config->getBinOffset(stage)] = true;
    }

    free(toSubtract);
}

int BackEnd::getHash(int stage, int spectralIndex)
{
    bool** subsamplingMatrix = config->getSubsamplingMatrix(stage);
    int ret = 0;

    for (int r = 0; r < config->getBinNumIndices(stage); r++) {
        bool currState = false;
        for (int c = 0; c < config->getSignalBitLength(); c++) {
            // currState = currState ^ (subsamplingMatrix[r][c] & (spectralIndexBinary[c] == 1));
            int shifted = 1 << (config->getSignalBitLength() - c - 1);
            currState = currState ^ (subsamplingMatrix[r][c] & ( (spectralIndex & shifted) == shifted));

        }
        if (currState)
            ret += 1<< (config->getBinNumIndices(stage)- r-1);
    }
    return ret;
}

void BackEnd::getClusteredFrequencies()
{
    // DecodedFrequencies is an unordered_map hence the elements are
    // not ordered with respect to its key
    // We order them here
    std::map<int,ffast_complex> decodedFrequenciesOrdered;
    
    double energy;
    double totalEnergy = 0;
    double peak = 0;
    ffast_complex amplitude=0;

    // Order the decoded frequencies with respect to their locations
    for (auto it = decodedFrequencies.cbegin(); it != decodedFrequencies.cend(); ++it)
    {
        decodedFrequenciesOrdered.emplace(it->first,it->second);
    }

    // Clear the decoded frequencies map
    decodedFrequencies.clear();

    // Ratio of the original signal length to the truncated one
    double ratio = ((double) config->getSignalLengthOriginal())/((double) config->getSignalLength());
    
    int nb = config->getSignalSparsity() * 100;

    for (auto it = decodedFrequenciesOrdered.cbegin(); it != decodedFrequenciesOrdered.cend(); ++it)
    {
        energy = std::norm(it->second);
        peak  += energy * (it->first);
        totalEnergy += energy;
        amplitude   += (it->second) * energy;
	// if the cluster has a band of two zeros afterwards
        if( decodedFrequenciesOrdered.count(it->first+1) + decodedFrequenciesOrdered.count(it->first+2) ==  0 )
        {
	    // weighted average location
            realFrequenciesIndices.push_back( ratio * peak / totalEnergy );
	    
	    // round to the nearest integer location
            decodedFrequencies[(int) round(ratio * peak / totalEnergy)] = std::polar(sqrt(totalEnergy),0.0);

            // Begin Average Method
            decodedFrequencies[(int) round(ratio * peak / totalEnergy)] = 0;
	    
            for (int i = 0; i < nb; ++i)
            {
                decodedFrequencies[(int) round(ratio*peak/totalEnergy)] += 
		    (frontEnd->getTimeSignal().at(i)) * std::polar(1.0,-2*M_PI*(i)*(peak/totalEnergy)/config->getSignalLength());
            }
            decodedFrequencies[(int) round(ratio*peak/totalEnergy)] = 
		std::polar(sqrt(totalEnergy),std::arg(decodedFrequencies[(int) round(ratio*peak/totalEnergy)]));
            // End Average Method

            amplitude = 0;
            totalEnergy = 0;
            peak = 0;
        }
    }
}
