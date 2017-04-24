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
    int iterations;
    int signalLength;
    int signalLengthOriginal;
    bool maximumLikelihood;
    int signalSparsity;
    std::vector<int> bins;
    std::vector<int> binOffsets;
    int binsSum;
    int delaysNb;
    int chainsNb;
    int delaysPerBunchNb;
    bool noisy;
    float SNRdB;
    int phasesNb;
    bool verbose;
    bool displayIterationTime;
    bool reconstructSignalInBackEnd;
    bool countSamples;
    std::vector<double> distribution;
    int signalSparsityPeeling;
    float maxSNRdB;
    bool defaultDelays;
    char* inputFile;
    char* outputFile;
    int signalLogLength;
    std::vector<bool**> subsamplingMatrices;

public:
    Config(int newArgc, char** newArgv);
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
    float getSNRdB() const;
    float getMaxSNRdB() const;
    int getPhasesNb() const;
    int getMaxErrorsBeforeQuitting() const;
    bool isVerbose() const; 
    bool needToDisplayIterationTime() const;
    bool needToReconstructSignalInBackEnd() const;
    bool needToCountSamples() const;
    void display() const;
    void help();
    bool isHelpDisplayed() const;
    std::vector<double> getDistribution() const;
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
};

#endif // CONFIG_H
