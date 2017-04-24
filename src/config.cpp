#include "config.h"
#include "split.h"
#include "utils.h"

#include <algorithm>
#include <assert.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <vector>



Config::Config(int newArgc, char** newArgv):argc(newArgc), argv(newArgv), helpDisplayed(false)
{
    bins.clear();
    setDefaultOptions();

    setOptionsFromCommandLine();

    if(distribution.size() == 0)
    {
      preprocessDistribution((char*) "1");
    }

    setBinOffsetsAndBinsSum();
    setSubsamplingMatrices();

    maxSNRdB = (maxSNRdB < SNRdB) ? SNRdB : maxSNRdB;


    chainsNb = (chainsNb >= 1) ? chainsNb : 1;
    if (!noisy)
      chainsNb = 1;

    delaysPerBunchNb = getSignalBitLength()+1;
    delaysNb = chainsNb * delaysPerBunchNb;
    


   assert( signalLength >= signalSparsity );
}

void Config::setDefaultOptions()
{
    outputFile = (char*) "sprightOutput.txt";
    
    signalLength = 64;
    signalLengthOriginal = signalLength;
    signalSparsityPeeling = 4;
    signalSparsity = 4;
    maximumLikelihood = false;
    countSamples = true;
    delaysPerBunchNb = getSignalBitLength() + 1;

    chainsNb = 1;

    /* for experiment mode */
    iterations = 1;
    experimentMode = true;
    // Phase = 0 implies phase of non-zero 
    // coefficients is uniformly random in [0,2*pi]. 
    phasesNb = 0; 
    displayIterationTime = false;
    noisy = false;
    SNRdB = 50;
    maxSNRdB = -std::numeric_limits<float>::infinity();
    verbose = false;
    reconstructSignalInBackEnd = false;
    defaultDelays = true;
    signalLogLength = log2(signalLength);

    bins.clear();
}

bool Config::isExperimentMode() const
{
    return experimentMode;
}

int Config::getIterations() const
{
    return iterations;
}

int Config::getSignalLength() const
{
    return signalLength;
}

bool Config::needToUseMaximumLikelihoodDetection() const
{
    return maximumLikelihood;
}

int Config::getSignalSparsity() const
{
    return signalSparsity;
}

int Config::getSignalSparsityPeeling() const
{
    return signalSparsityPeeling;
}

int Config::getSamplingPattern(int stage) const
{
    return bins[stage];
}

int Config::getBinNumItems(int stage) const
{
    return pow(2,numOnes(bins[stage]));
}

int Config::getBinNumIndices(int stage) const 
{
    return numOnes(bins[stage]);
}

int Config::getBinOffset(int stage) const
{
    return binOffsets[stage];
}

int Config::getBinsSum() const
{
    return binsSum;
}

int Config::getBinsNb() const
{
    return bins.size();
}

int Config::getDelaysNb() const
{
    return delaysNb;
}

int Config::getChainsNb() const
{
    return chainsNb;
}

int Config::getDelaysPerBunchNb() const
{
    return delaysPerBunchNb;
}

bool Config::isNoisy() const
{
    return noisy;
}

float Config::getSNRdB() const
{
    return SNRdB;
}

float Config::getMaxSNRdB() const
{
    return maxSNRdB;
}

int Config::getPhasesNb() const
{
    return phasesNb;
}


bool Config::isVerbose() const
{
    return verbose;
}

bool Config::needToDisplayIterationTime() const
{
    return displayIterationTime;
}

bool Config::needToReconstructSignalInBackEnd() const
{
    return reconstructSignalInBackEnd;
}

bool Config::needToCountSamples() const
{
    return countSamples;
}

int Config::getSignalBitLength() const
{
    return signalLogLength;
}

bool** Config::getSubsamplingMatrix(int stage) const
{
    return subsamplingMatrices[stage];
}

void Config::display() const
{
    std::cout << std::endl << "<=======================>" << std::endl;
    
    if( isExperimentMode() )
    {

       std::cout << "Running experiment mode" << std::endl;
   }

   std::cout << std::endl << "<===== INPUT PROFILE =====>" << std::endl;

   std::cout << signalLength << " -> signal length" << std::endl;
   std::cout << signalSparsity << " -> signal sparsity" << std::endl;
   if (noisy)
   {
    std::cout << SNRdB << " dB" << " -> signal-to-noise ratio (SNR)" << std::endl;
}
else
{
    std::cout << "Noiseless signal" << std::endl;
}
if (phasesNb == 0)
{
    std::cout << "Random phase"  << std::endl;
}
else
{
    std::cout << phasesNb << " -> possible phase";

    if (phasesNb > 1)
    {
        std::cout << "s";
    }

    std::cout << " for the non-zero frequencies in the input signal" << std::endl;
}

std::cout << std::endl << "<===== SPRIGHT CONFIG =====>" << std::endl;

std::cout << chainsNb << " -> delay chain(s)" << std::endl;

if( signalLengthOriginal != signalLength )
{
	std::cout << "The signal length is not a power of 2." << std::endl; 
	std::cout << "The last " << (signalLengthOriginal - signalLength) << " samples will be deleted" << std::endl;
}

for (unsigned int stage=0; stage<bins.size(); stage++)
{
    std::cout << bins[stage] << " ";
}
std::cout << "-> bins" << std::endl;

if (iterations > 1)
{
    std::cout << iterations << " -> iterations"  << std::endl;
}





if (maximumLikelihood)
{
    std::cout << "Maximum likelihood detection enabled"  << std::endl;
}

if (reconstructSignalInBackEnd)
{
    std::cout << "The back-end will reconstruct the full length frequency signal" << std::endl;
}

if (verbose)
{
    std::cout << "Verbose mode enabled" << std::endl;

    if (displayIterationTime)
    {
        std::cout << "Prompting execution time for each iteration" << std::endl;
    }
}

if(distribution.size() > 2) 
{
    std::cout << "Non uniform distribution -> ";

    for (unsigned int i = 0; i < distribution.size(); ++i)
    {
        std::cout << " " << distribution[i];
    }

    std::cout << std::endl;
}
}

void Config::help()
{
    helpDisplayed = true;

    std::cout   << std::endl << "<===== HELP =====>"
    << std::endl << std::endl
    << " [-a or --experiment]" << std::endl
    << "     Run experiment mode."  
    << std::endl << std::endl
    << " [-b BINS or --bins BINS]" << std::endl
    << "     Set the bins to use." << std::endl
    << "     Example: -b \"3 12 48\"."
    << std::endl << std::endl
    << " [-c or --samples]" << std::endl
    << "     Do not count the number of signal samples used in the front-end and display it in the results section."
    << std::endl << std::endl
    << " [-e NUM or --chains NUM]" << std::endl
    << "     Set the number of chains to use in the back-end for fast search."
    << std::endl << std::endl
    << " [-f FNAME or --file FNAME]" << std::endl
    << "     Input file to use" << std::endl
    << "     Example: -f \"timeSignal.txt\"."
    << std::endl << std::endl
    << " [-h or --help]" << std::endl
    << "     Displays help."
    << std::endl << std::endl
    << " [-i NUM or --iterations NUM]" << std::endl
    << "     Set the number of iterations." << std::endl
    << "     When the algorithm is be executed several times, the results will be averaged."
    << std::endl << std::endl
    << " [-k NUM or --sparsity NUM]" << std::endl
    << "     Set the signal sparsity. The input signal will have k uniformly distributed non zero frequency." << std::endl
    << "     Note: the following assertions should be true (k <= n) and (k <= 2*min(Bins))."
    << std::endl << std::endl
    << " [-l or --ml]" << std::endl
    << "     Enable the maximum likelihood detection." << std::endl
    << "     Note: This may be required when the snr is really low but it will dramatically slow down the execution."
    << std::endl << std::endl
    << " [-n NUM or --length NUM]" << std::endl
    << "     Set the signal length."
    << std::endl << std::endl
    << " [-r or --reconstruct]" << std::endl
    << "     Force the back-end to reconstruct the full n-length frequency signal." << std::endl
    << "     Note: this can increase the execution time of the SPRIGHT algorithm."
    << std::endl << std::endl
    << " [-s NUM or --snr NUM]" << std::endl
    << "     Add noise to the input signal in frequency domain and specify the SNR in dB." << std::endl
    << "     Example: -s 15 will generate a noisy signal with a SNR equal to 15 dB."
    << std::endl << std::endl
    << " [-v or --verbose]" << std::endl
    << "     Enable the verbose mode: the details of each iteration will be displayed."
    << std::endl << std::endl
    << " [-u DIST or --distribution DIST]" << std::endl
    << "     Enable the non-uniform distribution mode: input the weights of your distribution." << std::endl
    << "     Example: -u \"3 2 1\" will divide your frequency domain in 3, the probability to be in the first third will be 1/2, the next 1/3 and the last 1/6"
    << std::endl << std::endl
    << " [-z FNAME or --write FNAME]" << std::endl
    << "     Write the output on the choosen file"
    << std::endl;
}

bool Config::isHelpDisplayed() const
{
    return helpDisplayed;
}

void Config::setOptionsFromCommandLine()
{

    const struct option longOptions[] =
    {	
       {"experiment",	 no_argument,	    NULL, 'a'},
       {"bins",	 required_argument, NULL, 'b'},
       {"samples",      no_argument,       NULL, 'c'},
       {"chains",       required_argument, NULL, 'e'},
       {"file",         required_argument, NULL, 'f'},
       {"write",        required_argument, NULL, 'w'},
       {"help",         no_argument,       NULL, 'h'},
       {"iterations",   required_argument, NULL, 'i'},
       {"sparsity",     required_argument, NULL, 'k'},
       {"ml",           no_argument,       NULL, 'l'},
       {"length",       required_argument, NULL, 'n'},
       {"optimize",     no_argument,       NULL, 'o'},
       {"reconstruct",  no_argument,       NULL, 'r'},
       {"snr",          required_argument, NULL, 's'},
       {"distribution", required_argument, NULL, 'u'},
       {"verbose",      no_argument,       NULL, 'v'},
       {0, 0, 0, 0}
   };

   int option = 0;

   while(option != -1)
   {
        option = getopt_long(argc, argv, "acf:e:hi:k:b:ln:op:rs:tu:vw:z:", longOptions, NULL); //abfgjuxyz
        switch (option)
        {
           case 'a':
           experimentMode = true;
           break;

           case 'b':
           setBins(optarg);
           break;

           case 'n':
           signalLengthOriginal = atoi(optarg);
           signalLogLength = floor(log2(signalLengthOriginal));
           signalLength = pow(2,signalLogLength);
           break;

           case 'c':
           countSamples = false;
           break;

           case 'e':
           chainsNb = atoi(optarg);
           defaultDelays = false;
           break;

           case 'f':
           inputFile = optarg;
           experimentMode = false;
           break;

           case 'h':
           help();
           break;

           case 'i':
           iterations = atoi(optarg);
           break;

           case 'k':
           signalSparsity = atoi(optarg);
           signalSparsityPeeling = atoi(optarg);
           break;

           case 'l':
           maximumLikelihood = true;
           break;

           case 'r':
           reconstructSignalInBackEnd = true;
           break;

           case 's':
           noisy = true;
           SNRdB = strtof(optarg, 0);
           break;

           case 'u':
           preprocessDistribution(optarg);
           break;

           case 'z':
           outputFile = optarg;
           break;

           case 'v':
           verbose = true;
           break;
       }
   }
}

void Config::setBins(const char *newBins)
{
    bins.clear();
    // split: divides a string into multiple strings that are separated by space ' '.
    std::vector<std::string> tempBins = split(newBins, ' ');
    bins.clear();
    for(unsigned int stage=0; stage<tempBins.size(); stage++)
    {
        int value = atoi(tempBins[stage].c_str());
        bins.push_back(value);
    }
}

void Config::setBinOffsetsAndBinsSum()
{
    if (bins.size() == 0) {
        bins.clear();
        int floorThird = (int) (signalLogLength/3);
        bins.push_back(pow(2,floorThird)-1);
        if (signalLogLength % 3 == 0) {
            bins.push_back(pow(2,2*floorThird)-pow(2,floorThird));
            bins.push_back(pow(2,3*floorThird)-pow(2,2*floorThird));
        }
        if (signalLogLength % 3 == 1) {
            bins.push_back(pow(2,2*floorThird)-pow(2,floorThird));
            bins.push_back(pow(2,3*floorThird+1)-pow(2,2*floorThird));
        }
        if (signalLogLength % 3 == 2) {
            bins.push_back(pow(2,2*floorThird+1)-pow(2,floorThird));
            bins.push_back(pow(2,3*floorThird+2)-pow(2,2*floorThird+1));    
        }
    }


    binsSum = 0;
    for (unsigned int stage=0; stage<bins.size(); stage++)
    {
        binOffsets.push_back(binsSum);
        binsSum += pow(2,numOnes(bins[stage]));
    }
}

void Config::setSubsamplingMatrices()
{
    subsamplingMatrices.clear();
    int samplingPattern;
    int r;
    for (int stage = 0; stage < getBinsNb(); stage++) {
        samplingPattern = getSamplingPattern(stage);

        bool** subsamplingMatrix = (bool**) malloc(getBinNumIndices(stage) * sizeof(bool*));

        for (int iter=0; iter<getBinNumIndices(stage); iter++)
        {
            subsamplingMatrix[iter] = (bool*) malloc(getSignalBitLength() * sizeof(bool));
        }
        for (int r = 0; r < getBinNumIndices(stage); r++)
            for (int c = 0; c < getSignalBitLength(); c++)
                subsamplingMatrix[r][c] = 0;

        r = 0;
        for (int c = 0; c < getSignalBitLength() & r < getBinNumIndices(stage); c++) {
            if ((samplingPattern & (1 << (getSignalBitLength() - c - 1))) != 0) {

                subsamplingMatrix[r++][c] = 1;
            }
        }
        subsamplingMatrices.push_back(subsamplingMatrix);
    }
}

std::vector<double> Config::getDistribution() const
{
  return distribution;
}

void Config::preprocessDistribution(char* newDistribution)
{
    std::vector<std::string> tempDistribution = split(newDistribution, ' ');

    distribution.clear();
    distribution.push_back(0);

    int l = tempDistribution.size()+1;

    for (int i = 1; i < l; ++i)
    {
       distribution.push_back(distribution[i-1]+atof(tempDistribution[i-1].c_str()));
   }

   for (int i = 1; i < l; ++i) 
   {
       distribution[i] = distribution[i]/distribution[l-1];
   }
}

int Config::getSignalLengthOriginal() const
{
    return signalLengthOriginal;
}

char* Config::getInputFile() const
{
    return inputFile;
}

char* Config::getOutputFile() const
{
    return outputFile;
}
