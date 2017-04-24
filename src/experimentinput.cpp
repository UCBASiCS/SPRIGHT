#include "input.h"
#include "experimentinput.h"
#include "utils.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>

ExperimentInput::ExperimentInput(Chrono* newChrono, const Config* newConfig): Input(newChrono, newConfig)
{
    // set the signal amplitudes to 1 and change
    // the noise variance to account for SNR
    signalMagnitude = 1;
    noiseStdDeviation = pow(10,-config->getSNRdB()/20);
}

ExperimentInput::~ExperimentInput()
{
}

void ExperimentInput::process()
{
    chrono->start("Input");
    neededSamples.clear();
    for (int k = 0; k < config->getSignalLength(); k++)
    {
	   neededSamples.insert(k);
    }

    nonZeroFrequencies.clear();
    generateNonZeroFrequencies();

    for (auto t = neededSamples.cbegin(); t != neededSamples.cend(); ++t) {
        timeSignal[*t] =0;
    }
    
    frequencyToTime();
    
    if (config->isNoisy())
    {
        addNoise();
    }
    
    for (auto t = neededSamples.cbegin(); t != neededSamples.cend(); ++t)
    {
        timeSignal[*t] /= sqrt((ffast_complex)config->getSignalLengthOriginal());
    }
    chrono->stop("Input");
}

void ExperimentInput::process(std::vector<int> delays)
{

    chrono->start("Input");
    findNeededSamples(delays);

    generateNonZeroFrequencies();
    

    // frequencyToTimeUsingFFT(delays);
    frequencyToTime();
        

    if (config->isNoisy())
    {
        addNoise();
    }

    for (auto t = neededSamples.cbegin(); t != neededSamples.cend(); ++t)
    {
        timeSignal[*t] /= sqrt((ffast_complex)config->getSignalLengthOriginal());
    }


    
    chrono->stop("Input");
}

const std::unordered_map<int,ffast_complex>& ExperimentInput::getTimeSignal() const
{
    return timeSignal;
}

ffast_complex ExperimentInput::getSignalAtIndex(int k) const 
{
    return timeSignal.at(k);
}

const std::unordered_map<int, ffast_complex> &ExperimentInput::getNonZeroFrequencies() const
{
    return nonZeroFrequencies;
}

ffast_real ExperimentInput::getSignalMagnitude() const
{
    return signalMagnitude;
}

ffast_real ExperimentInput::getRealSNR() const
{
    return realSNR;
}

ffast_real ExperimentInput::getRealSNRdB() const
{
    return (ffast_real) 10*log10((double) realSNR);
}

ffast_real ExperimentInput::distribution(ffast_real _urand, std::vector<double> F) 
{
    int l = F.size();
    for (int i = 1; i < l; ++i)
    {
        if(_urand < F[i])
            return ( ( (_urand-F[i])/(F[i]-F[i-1]) + ( (double) i) )/( (double) l-1.0) );
    }
    return 1;
}

void ExperimentInput::generateNonZeroFrequencies()
{
    std::set<int> tempLocations;
    nonZeroFrequencies.clear();

    while((int) tempLocations.size() < config->getSignalSparsity())
    {
        int tempLocation = ((int) floor(config->getSignalLengthOriginal()*distribution((ffast_real) drand48(), config->getDistribution())))
                % config->getSignalLengthOriginal();
        tempLocations.insert(tempLocation);
    }

    for(auto it = tempLocations.cbegin(); it != tempLocations.cend(); ++it)
    {
        nonZeroFrequencies[*it] = std::polar(signalMagnitude, getRandomPhase());
    }
}

void ExperimentInput::printNonZeroFrequencies() const
{
    for ( auto it = nonZeroFrequencies.begin(); it != nonZeroFrequencies.end(); ++it )
        std::cout << it->first << ":" << it->second <<", ";
    std::cout << std::endl;
}

ffast_real ExperimentInput::getRandomPhase() const
{
    if (config->getPhasesNb() < 1)
    {
       return 2*M_PI*drand48();
    }
    else if (config->getPhasesNb() == 1)
    {
        return 0;
    }
    else if (config->getPhasesNb() == 2)
    {
        return M_PI*std::round(drand48());
    }
    else
    {
        return (std::floor(config->getPhasesNb()*drand48())*2 + 1)*M_PI/config->getPhasesNb();
    }
}

void ExperimentInput::addNoise()
{
    ffast_real signalPower = signalMagnitude * signalMagnitude;
    ffast_real noisePower = 0;
    ffast_real noiseNormFactor;
    ffast_real noisePhase;

    for (auto t = neededSamples.cbegin(); t != neededSamples.cend(); ++t)
    {
        noiseNormFactor = (ffast_real) -2*log(drand48());
        noisePhase = (ffast_real) 2*M_PI*drand48();
        noisePower += noiseNormFactor; // Standard deviation added later.
        timeSignal[*t] += std::polar( noiseStdDeviation * ((ffast_real) sqrt( (double) noiseNormFactor ) ) / sqrt(2) , noisePhase );
    }

    noisePower *= noiseStdDeviation * noiseStdDeviation;
    realSNR = signalPower * config->getSignalLength()/noisePower;
}

void ExperimentInput::frequencyToTime()
{
    // go through the required sample indices
    for (auto t = neededSamples.cbegin(); t != neededSamples.cend(); ++t)
    {
        timeSignal[*t] = 0;
    	// go through the non-zero frequency components
        for (auto f = nonZeroFrequencies.cbegin(); f != nonZeroFrequencies.cend(); ++f)
        {
            timeSignal[*t] += f->second * (numOnes(*t & f->first) % 2 ? -1.0 : 1.0);
        }
    }
}

void ExperimentInput::frequencyToTimeUsingFWHT(std::vector<int> delays)
{
    // unimplemented
}

void ExperimentInput::findNeededSamples(std::vector<int> delays)
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
