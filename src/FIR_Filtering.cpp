#include "FIR_Filtering.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

FIR_Filtering::FIR_Filtering()
{
}

bool FIR_Filtering::getWindowedSincFIRLPCoefficients(double dSamplingFrequency, double dCutOffFrequency, unsigned int nFilterOrder, std::vector<double> &firCoefficients)
{
	// check valid filter parameters
	if (dSamplingFrequency <= 0.0 || dCutOffFrequency <= 0.0)
		return false;

	if (dCutOffFrequency / dSamplingFrequency > 0.5)
		return false;

	// relative cutoff frequency
	double FC = dCutOffFrequency / dSamplingFrequency;

	// filter order
	unsigned int M = nFilterOrder;

	// minimum order M = 2
	if (M < 2)
		return false;

	firCoefficients.resize(M);

	for (unsigned int i = 0; i < M; i++)
	{
		// generate si = sin(x)/x function values
		if (i - M / 2 == 0)
		{
			firCoefficients[i] = 2.0*M_PI*FC;
		}
		else
		{
			firCoefficients[i] = sin(2.0*M_PI*FC * (double(i) - double(M) / 2.0)) / (double(i) - double(M) / 2.0);
		}

		// apply hanning window
		firCoefficients[i] = firCoefficients[i] * (0.5* (1.0 - cos(2.0*M_PI*(double(i)) / double(M - 1))));
	}

	// normalize filter coefficients
	double dSum = 0.0;
	for (unsigned int i = 0; i < M; i++)
	{
		dSum += firCoefficients[i];
	}

	if (dSum == 0.0)
		return false;

	for (unsigned int i = 0; i < M; i++)
	{
		firCoefficients[i] = firCoefficients[i] / dSum;
	}

	return true;
}

bool FIR_Filtering::getWindowedSincFIRHPCoefficients(double dSamplingFrequency, double dCutOffFrequency, unsigned int nFilterOrder, std::vector<double> &firCoefficients)
{
	// check valid filter parameters
	if (dSamplingFrequency <= 0.0 || dCutOffFrequency <= 0.0)
		return false;

	if (dCutOffFrequency / dSamplingFrequency > 0.5)
		return false;

	// get first low pass filter coefficients
	if (!getWindowedSincFIRLPCoefficients(dSamplingFrequency, dCutOffFrequency, nFilterOrder, firCoefficients))
		return false;

	// filter order
	unsigned int M = nFilterOrder;

	// change the low pass filter kernel to a high pass filter kernel by
	// spectral inverision
	for (unsigned int i = 0; i < M; i++)
	{
		firCoefficients[i] = -firCoefficients[i];
	}

	// add one to coefficient in center of filter kernel
	firCoefficients[M / 2] = firCoefficients[M / 2] + 1.0;

	return true;
}

// static function, used to perform filtering without instantiation of filter class
bool FIR_Filtering::filter1D(const std::vector<double> &inputSignal, std::vector<double> &outputSignal, std::vector<double> firCoefficients)
{

	// allocate output vector
	outputSignal.resize(inputSignal.size());

#pragma omp parallel for// private (data, kernel)
	for (long lSample = 0; lSample < outputSignal.size(); lSample++)
	{
		long nFilterLength = firCoefficients.size();
		if (nFilterLength > lSample)
			nFilterLength = lSample + 1;

		for (long k = 0; k < nFilterLength; k++)
		{
			outputSignal[lSample] += firCoefficients[k] * inputSignal[lSample - k];
		}

		// adapt first samples
		for (long k = nFilterLength; k < firCoefficients.size(); k++)
		{
			outputSignal[lSample] += firCoefficients[k] * inputSignal[0];
		}
	}

	return true;
}