#include "input.h"
#include "customizedinput.h"
#include "utils.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>

CustomizedInput::CustomizedInput(Chrono* newChrono, const Config* newConfig): Input(newChrono, newConfig)
{}

CustomizedInput::~CustomizedInput()
{}

void CustomizedInput::process()
{
    chrono->start("Input");

    for (int k = 0; k < config->getSignalLength(); k++)
    {
	neededSamples.insert(k);
    }

    convertFileToTimeSignal();
    
    chrono->stop("Input");
}

void CustomizedInput::process(std::vector<int> delays)
{
    chrono->start("Input");

    findNeededSamples(delays);
    convertFileToTimeSignal();

    chrono->stop("Input");
}

const std::unordered_map<int,ffast_complex>& CustomizedInput::getTimeSignal() const
{
    return timeSignal;
}

ffast_complex CustomizedInput::getSignalAtIndex(int k) const 
{
    return timeSignal.at(k);
}


void CustomizedInput::convertFileToTimeSignal()
{
    std::string line;
    std::ifstream file(config->getInputFile());
    ffast_complex temp;

    for(int i=0; i < config->getSignalLength(); i++) {
        if(neededSamples.count(i) > 0 && getline(file,line)) {
            (std::istringstream) line >> temp;
            timeSignal[i] = temp;
        }
        else if (neededSamples.count(i) == 0 && file.ignore(256,'\n'))
        {}
        else {
            std::cerr << "Error while reading input" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

void CustomizedInput::findNeededSamples(std::vector<int> delays)
{

    int signalBitLength = config->getSignalBitLength();
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
                int mlInt = mapToInt(bin_iter, b, binIndices, signalBitLength);
                int signalIndex = currDelay ^ mlInt;

                neededSamples.insert(signalIndex);
            }
        }     
    }
}

void CustomizedInput::printNonZeroFrequencies() const
{
}