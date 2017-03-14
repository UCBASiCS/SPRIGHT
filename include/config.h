#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

class Config
{
private:
    int argc;
    char** argv;

    bool experimentMode;
    bool helpDisplayed;
    bool offGrid;
    int iterations;
    int signalLength;
    int signalLengthOriginal;
    bool maximumLikelihood;
    int lengthFactor;
    int signalSparsity;
    std::vector<int> bins;
    std::vector<int> binOffsets;
    int binsSum;
    int delaysNb;           // in draft paper : N=CD
    int chainsNb;           // in draft paper : C=noisy ? log(log(n)) : 1
    int delaysPerBunchNb;   // in draft paper : D=input[def:2]
    bool noisy;
    float SNRdB;
    float minFourierMagnitude;
    float offgridSNRdB;
    int phasesNb;
    bool verbose;
    bool compareWithFFTW;
    bool displayIterationTime;
    bool reconstructSignalInBackEnd;
    bool countSamples;
    int FFTWstrategy;
    std::vector<double> distribution;
    bool applyWindowVar;
    int signalSparsityPeeling;
    float maxSNRdB;
    bool defaultDelays;
    char* inputFile;
    char* outputFile;
    int signalLogLength;
    std::vector<bool**> subsamplingMatrices;
    // std::vector<int> randomVectorSampleName;

    // int hello;

public:
    Config(int newArgc, char** newArgv);
    bool isOffTheGrid() const;
    bool isExperimentMode() const;
    int getIterations() const;
    int getSignalLength() const;
    bool needToUseMaximumLikelihoodDetection() const;
    int getSignalSparsity() const;
    int getSignalSparsityPeeling() const;
    int getSamplingPattern(int stage) const;
    int getBinNumItems(int stage) const;
    int getBinNumIndices(int stage) const;
    int getBinOffset(int stage) const;
    int getBinsSum() const;
    int getBinsNb() const;
    int getDelaysNb() const;
    int getChainsNb() const;
    int getDelaysPerBunchNb() const;
    bool isNoisy() const;
    float getOffgridSNRdB() const;
    float getMinFourierMagnitude() const;
    float getSNRdB() const;
    float getMaxSNRdB() const;
    int getPhasesNb() const;
    int getMaxErrorsBeforeQuitting() const;
    bool isVerbose() const; 
    bool needToCompareWithFFTW() const;
    bool needToDisplayIterationTime() const;
    bool needToReconstructSignalInBackEnd() const;
    bool needToCountSamples() const;
    int getFFTWstrategy() const;
    void display() const;
    void help();
    bool isHelpDisplayed() const;
    std::vector<double> getDistribution() const;
    bool applyWindow() const;
    int getSignalLengthOriginal() const;
    char* getInputFile() const;
    char* getOutputFile() const;
    int getSignalBitLength() const;
    bool** getSubsamplingMatrix(int stage) const;


private:
    void setDefaultOptions();
    void setOptionsFromCommandLine();
    void setBinOffsetsAndBinsSum();
    void setSubsamplingMatrices();
    void setBins(const char* newBins);
    void setMaxErrorsBeforeQuitting();
    void preprocessDistribution(char* newDistribution);
    void proposedBins(int newRangeSemiLength = 10);
};

#endif // CONFIG_H
