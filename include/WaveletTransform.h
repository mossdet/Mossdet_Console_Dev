#include <vector>

#include "DataTypes.h"

#pragma once
class WaveletTransform
{
public:
	WaveletTransform(double samplingRate, double startFreq, double endFreq, double deltaFreq, int nrOscills);
	~WaveletTransform();

	bool run(matrixStd &data);
	unsigned int getNumberComponents() const;
	std::vector<double> getConvolutionKernel(double dCenterFrequency);
	bool getResult(matrixStd &outParameters);

protected:
	double m_dSamplingRate;
	double m_dStartFrequency;
	double m_dEndFrequency;
	double m_dDeltaFrequency;
	int m_nrOscillations;
	matrixStd m_Parameters;
};

