#ifndef SPRIGHT_H
#define SPRIGHT_H

#include "backend.h"
#include "chrono.h"
#include "config.h"
#include "fftw.h"
#include "frontend.h"
#include "input.h"
#include "output.h"

class SPRIGHT
{
public:
    const Config* config;

private:
    Chrono* chrono;
    Input* input;
    Output* output;

    int iteration;

public:
    SPRIGHT(const Config* newConfig, Input* newInput, Output* newOutput);
    ~SPRIGHT();
    void process();
    void displayResults() const;
    const std::vector<int> getDelays() const;
    BackEnd* backEnd;
    FrontEnd* frontEnd;

private:
    void displayExecutionTime(double frontEndTime, double backEndTime, double globalTime) const;
};

#endif // SPRIGHT_H
