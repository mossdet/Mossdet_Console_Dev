#pragma once

#include <vector>

class FIR_Filtering
{
public:
	FIR_Filtering();
	static bool getWindowedSincFIRLPCoefficients(double dSamplingFrequency, double dCutOffFrequency, unsigned int nFilterOrder, std::vector<double> &firCoefficients);
	static bool getWindowedSincFIRHPCoefficients(double dSamplingFrequency, double dCutOffFrequency, unsigned int nFilterOrder, std::vector<double> &firCoefficients);
	static bool filter1D(const std::vector<double> &inputSignal, std::vector<double> &outputSignal, std::vector<double> firCoefficients);
};

