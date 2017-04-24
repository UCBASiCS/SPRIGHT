#include "frontend.h"
#include "utils.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>

FrontEnd::FrontEnd(Chrono* newChrono, const Config* newConfig, const Input* newInput):
	Step(newChrono, newConfig),countSamplesDone(false),input(newInput)
{
    signalLength = config->getSignalLength();
    observationMatrix = (ffast_complex**) malloc(config->getBinsSum() * sizeof(ffast_complex*));

    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
    {
        observationMatrix[binAbsoluteIndex] = (ffast_complex*) malloc(config->getDelaysNb() * sizeof(ffast_complex));
    }
    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
        for (int delayNb = 0; delayNb < config->getBinsSum(); delayNb++)
            observationMatrix[binAbsoluteIndex][delayNb] = 0;

    computeDelays();
}

FrontEnd::~FrontEnd()
{
    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
    {
        free(observationMatrix[binAbsoluteIndex]);
    }

    free(observationMatrix);
}

void FrontEnd::process()
{
    chrono->start("FrontEnd");
    
    int realSignalBitLength = config->getSignalBitLength(); //ceil(log2(config->getSignalLength()))

    for (int binAbsoluteIndex=0; binAbsoluteIndex<config->getBinsSum(); binAbsoluteIndex++)
        for (int delayNb = 0; delayNb < config->getDelaysNb(); delayNb++)
            observationMatrix[binAbsoluteIndex][delayNb] = 0;

    for (int c = 0; c < config->getBinsNb(); c++)
    {
        int binIndices = config->getSamplingPattern(c);

        int b = config->getBinNumIndices(c);
        int B = pow(2,b); // number of bins in stage

        for (int delayIndex = 0; delayIndex < config->getDelaysNb(); delayIndex++)
        {
            int currDelay = delays[delayIndex];

            std::vector<ffast_complex> binTimeDomain = std::vector<ffast_complex>(B);
            for (int bin_iter = 0; bin_iter < B; bin_iter++)
            {
                int mlInt = mapToInt(bin_iter, b, binIndices, realSignalBitLength);
                int signalIndex = currDelay ^ mlInt;

                binTimeDomain[bin_iter] = input->getSignalAtIndex(signalIndex);

                if (!countSamplesDone && config->needToCountSamples())
                    usedSamples.insert(signalIndex);
            }


            std::vector<ffast_complex> a = std::vector<ffast_complex>(binTimeDomain);
            std::vector<ffast_complex> b = std::vector<ffast_complex>(B);
            std::vector<ffast_complex> tmp = std::vector<ffast_complex>(B);
            int ii, jj, ss;
            for (ii = B>>1; ii > 0; ii>>=1) {
                for (jj = 0; jj < B; jj++) {
                    ss = jj/ii%2;
                    b[jj]=a[(ss?-ii:0)+jj];
                    b[jj]+= ((ffast_complex) (ss?-1:1)) * a[(ss?0:ii)+jj];
                }
                tmp = a; a = b; b = tmp;
            }
            std::vector<ffast_complex> binTransform = std::vector<ffast_complex>(a);

            for (int binLoc = 0; binLoc < B; binLoc++)
                observationMatrix[config->getBinOffset(c)+binLoc][delayIndex] = binTransform[binLoc];
        }     
    }






    for (int stage = 0; stage < config->getBinsNb(); stage++) {
        for (int r = 0; r < config->getBinNumItems(stage); r++) {
            for (int c = 0; c < config->getDelaysNb(); c++) {
                ffast_complex val = observationMatrix[r + config->getBinOffset(stage)][c];
                // double multiplier = pow(1.0*config->getSignalLength() / config->getBinNumItems(stage),.5);
                double multiplier = 1.0/config->getBinNumItems(stage) * pow(config->getSignalLength(), .5);
                // double multiplier = 1.0/config->getBinNumItems(stage);

                observationMatrix[r + config->getBinOffset(stage)][c] = val * multiplier;
            }
        }
    }

    countSamplesDone = true;

    chrono->stop("FrontEnd");
}

ffast_complex** FrontEnd::getObservationMatrix() const
{
    return observationMatrix;
}

int FrontEnd::getUsedSamplesNb() const
{
    return usedSamples.size();
}

std::vector<int> FrontEnd::getDelays() const
{
    return delays;
}

void FrontEnd::computeDelays()
{
    delays.clear();

    if (config->isNoisy()) {
        if (config->needToUseMaximumLikelihoodDetection()) {
            if (config->isVerbose())
                std::cout << "ML delays" << std::endl;
            std::set<int> tempDelaysSet;
            while((int) tempDelaysSet.size() < config->getDelaysNb())
            {
                int tempDelay = ((int) floor(signalLength*((ffast_real) drand48()))) % signalLength;

                tempDelaysSet.insert(tempDelay);
            }
            delays.resize(config->getDelaysNb());
            std::copy(tempDelaysSet.begin(), tempDelaysSet.end(), delays.begin());
        }
        else {
            if (config->isVerbose())
                std::cout<< "FFAST chain delays" << std::endl;
            std::set<int> tempDelaysSet;
            while((int) tempDelaysSet.size() < config->getChainsNb())
            {
                int tempDelay = ((int) floor(signalLength*((ffast_real) drand48()))) % signalLength;

                tempDelaysSet.insert(tempDelay);
            }

            std::vector<int> tempDelaysArr;
            tempDelaysArr.clear();
            tempDelaysArr.resize(config->getDelaysNb());
            std::copy(tempDelaysSet.begin(), tempDelaysSet.end(), tempDelaysArr.begin());

            int r = 0;
            for (int i = 0; i < config->getChainsNb(); ++i) {
                r = tempDelaysArr[i];
                delays.push_back(r);
                for (int j = 0; j < config->getDelaysPerBunchNb()-1; j++)
                    delays.push_back(positiveMod(r^(1<<(config->getSignalBitLength() -1 - j)), signalLength));
            }
        }
    }
    else 
    {
        if (config->isVerbose())
            std::cout << "using standard delays" << std::endl;

        delays.push_back(0);
        for (int i =0; i < config->getSignalBitLength(); i++) {
            delays.push_back(1<<(config->getSignalBitLength() -1 - i));
        }
    }
}

const std::unordered_map<int,ffast_complex>& FrontEnd::getTimeSignal() const
{
    return input->getTimeSignal();
}
