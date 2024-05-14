#include <fstream>
#include <math.h>
#include <algorithm>
#include <direct.h>
#include <string>

#include <FIR_Filtering.h>
#include "WaveletTransform.h"
#include "SignalProcessing.h"
#include "Plotter.h"

#include <omp.h>

SignalProcessing::SignalProcessing()
{
}


SignalProcessing::~SignalProcessing()
{
}

bool SignalProcessing::spindleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate) {

	//Filter Signal with "samplesToRead" number of samples and create new matrixStd to evaluate
	std::vector <double> firCoefficients;
	int nrRows = signalSamples.size();
	int nrColumns = signalSamples[0].size();
	unsigned int M = FIR_RIPPLE_ORDER; // filter order


	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dHighCutOffFreq, M, firCoefficients);
	//Allocate
	matrixStd matrixA;
	matrixA.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		matrixA[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			matrixA[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}


	firCoefficients.clear();
	FIR_Filtering::getWindowedSincFIRHPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	//Allocate
	filteredSignalSamplesBP.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		filteredSignalSamplesBP[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(matrixA[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			filteredSignalSamplesBP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}

	return true;
}

bool SignalProcessing::rippleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate) {

	//Filter Signal with "samplesToRead" number of samples and create new matrixStd to evaluate
	std::vector <double> firCoefficients;
	int nrRows = signalSamples.size();
	int nrColumns = signalSamples[0].size();
	unsigned int M = FIR_RIPPLE_ORDER; // filter order


	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dHighCutOffFreq, M, firCoefficients);
	//Allocate
	matrixStd matrixA;
	matrixA.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		matrixA[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			matrixA[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}


	firCoefficients.clear();
	FIR_Filtering::getWindowedSincFIRHPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	//Allocate
	filteredSignalSamplesBP.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		filteredSignalSamplesBP[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(matrixA[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			filteredSignalSamplesBP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}

	return true;
}

bool SignalProcessing::fastRippleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate) {

	//Filter Signal with "samplesToRead" number of samples and create new matrixStd to evaluate
	std::vector <double> firCoefficients;
	int nrRows = signalSamples.size();
	int nrColumns = signalSamples[0].size();
	unsigned int M = FIR_FR_ORDER; // filter order


	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dHighCutOffFreq, M, firCoefficients);
	//Allocate
	matrixStd matrixA;
	matrixA.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		matrixA[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			matrixA[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}


	firCoefficients.clear();
	FIR_Filtering::getWindowedSincFIRHPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	//Allocate
	filteredSignalSamplesBP.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		filteredSignalSamplesBP[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(matrixA[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			filteredSignalSamplesBP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}

	return true;
}
bool SignalProcessing::applyBandPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate) {
	
	//Filter Signal with "samplesToRead" number of samples and create new matrixStd to evaluate
	std::vector <double> firCoefficients;
	int nrRows = signalSamples.size();
	int nrColumns = signalSamples[0].size();
	unsigned int M = FIR_ORDER; // filter order


	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dHighCutOffFreq, M, firCoefficients);
	//Allocate
	matrixStd matrixA;
	matrixA.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		matrixA[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			matrixA[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}


	firCoefficients.clear();
	FIR_Filtering::getWindowedSincFIRHPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	//Allocate
	filteredSignalSamplesBP.resize(nrRows);
	for (int i = 0; i < nrRows; ++i) {
		filteredSignalSamplesBP[i].resize(nrColumns);
	}
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrColumns);
		FIR_Filtering::filter1D(matrixA[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrColumns; sampleIndex++) {
			filteredSignalSamplesBP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrColumns - 1 - sampleIndex];
		}
	}

	return true;
}

bool SignalProcessing::applyHighPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesHP, double dLowCutOffFreq, double samplingRate) {

	int nrRows = signalSamples.size();
	int nrSamples;
	unsigned int M = FIR_HP_ORDER; // filter order
	std::vector <double> firCoefficients;
	//Allocate
	filteredSignalSamplesHP.resize(nrRows);
	for (int ch = 0; ch < nrRows; ++ch) {
	filteredSignalSamplesHP[ch].resize(signalSamples[ch].size());
	nrSamples = filteredSignalSamplesHP[ch].size();
	}

	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrSamples);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			filteredSignalSamplesHP[channIdx][sampleIndex] = signalSamples[channIdx][sampleIndex] - oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
	}

	/*int nrRows = signalSamples.size();
	int nrSamples;
	unsigned int M = FIR_ORDER*20; // filter order
	std::vector <double> firCoefficients;
	//Allocate
	filteredSignalSamplesHP.resize(nrRows);
	for (int ch = 0; ch < nrRows; ++ch) {
		filteredSignalSamplesHP[ch].resize(signalSamples[ch].size());
		nrSamples = filteredSignalSamplesHP[ch].size();
	}

	FIR_Filtering::getWindowedSincFIRHPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);
	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrSamples);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			filteredSignalSamplesHP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
	}*/

	return true;
}

bool SignalProcessing::applyLowPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesHP, double dLowCutOffFreq, double samplingRate) {

	int nrRows = signalSamples.size();
	int nrSamples;

	unsigned int M = FIR_ORDER; // filter order

	std::vector <double> firCoefficients;

	//Allocate
	filteredSignalSamplesHP.resize(nrRows);
	for (int ch = 0; ch < nrRows; ++ch) {
		filteredSignalSamplesHP[ch].resize(signalSamples[ch].size());
		nrSamples = filteredSignalSamplesHP[ch].size();
	}


	FIR_Filtering::getWindowedSincFIRLPCoefficients(samplingRate, dLowCutOffFreq, M, firCoefficients);

	for (int channIdx = 0; channIdx < nrRows; channIdx++) {
		std::vector<double> oneChannelFilteredSignal;
		std::vector<double> reversedSignal(nrSamples);
		FIR_Filtering::filter1D(signalSamples[channIdx], oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			reversedSignal[sampleIndex] = oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
		oneChannelFilteredSignal.clear();
		FIR_Filtering::filter1D(reversedSignal, oneChannelFilteredSignal, firCoefficients);
		for (long sampleIndex = 0; sampleIndex < nrSamples; sampleIndex++) {
			filteredSignalSamplesHP[channIdx][sampleIndex] = oneChannelFilteredSignal[nrSamples - 1 - sampleIndex];
		}
	}

	return true;
}

bool SignalProcessing::applyBandStopFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBS, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate) {
	int nrRows = signalSamples.size();
	int nrColumns;

	//Allocate
	filteredSignalSamplesBS.resize(nrRows);
	for (int ch = 0; ch < nrRows; ++ch) {
		filteredSignalSamplesBS[ch].resize(signalSamples[ch].size());
		nrColumns = filteredSignalSamplesBS[ch].size();
	}


	matrixStd stopSignal;
	SignalProcessing::applyBandPassFilterToSignal(signalSamples, stopSignal, dLowCutOffFreq, dHighCutOffFreq, samplingRate);

	for (unsigned r = 0; r < nrRows; ++r) {
		for (long c = 0; c < nrColumns; ++c) {
			filteredSignalSamplesBS[r][c] = signalSamples[r][c] - stopSignal[r][c];
		}
	}

	return true;
}

void SignalProcessing::getMatrixPower(matrixStd& matrixStdOutput) {
	int nrRows = matrixStdOutput.size();

	for (int r = 0; r < nrRows; r++) {
		long nrColumns = matrixStdOutput[r].size();
		for (long c = 0; c < nrColumns; c++) {
			matrixStdOutput[r][c] = matrixStdOutput[r][c] * matrixStdOutput[r][c];
		}
	}
}

bool SignalProcessing::firstDerivative(std::vector<double>& signal, std::vector<double>& firstDerivative, double startSample, double sampleLength, double samplingRate) {
	double step = 1.0 / samplingRate;
	firstDerivative.resize(sampleLength);
//#pragma omp parallel for
	for (long i = startSample; i < startSample + sampleLength; ++i) {
		if (i == 0) {
			firstDerivative[i - startSample] = (signal[i + 1] - signal[i]) / (step);
		}
		else if (i + 1 >= signal.size()) {
			firstDerivative[i - startSample] = (signal[i] - signal[i - 1]) / (step);
		}
		else {
			firstDerivative[i - startSample] = (signal[i + 1] - signal[i - 1]) / (2 * step);
		}
	}
	return true;
}


unsigned SignalProcessing::morletWaveletTransform(matrixStd &inData, double dSamplingFreq, double dStartFreq, double dEndFreq, double dDeltaFreq, int waveletOscills, matrixStd& waveletOutput) {

	// initial check vor validity
	if (dSamplingFreq< 1e-12)
		return false;

	unsigned int nNumberComponents = 0;
	WaveletTransform  wavelet(dSamplingFreq, dStartFreq, dEndFreq, dDeltaFreq, waveletOscills);
	wavelet.run(inData);
	// get result
	wavelet.getResult(waveletOutput);

	// get number components
	nNumberComponents = wavelet.getNumberComponents();

	return nNumberComponents;
}

bool SignalProcessing::setVarianceData(const std::vector<double> data) {
	if (data.size() == 0)
		return false;
	m_varData = data;
	return m_varData.size() == data.size();
}

double SignalProcessing::runVariance(long startSample, long sampleLength) {
	double variance = 0.0, average = 0.0;

	if (sampleLength < 2)
		return 0;

	for (int i = startSample; i < startSample + sampleLength; i++) {
		average += m_varData.at(i);
	}
	average /= sampleLength;

	for (int i = startSample; i < startSample + sampleLength; i++) {
		variance += ((m_varData.at(i) - average)*(m_varData.at(i) - average));
	}

	variance /= sampleLength - 1;

	return variance;
}

double SignalProcessing::getPearsonCorrCoeff(std::vector<double>& signalA, std::vector<double>& signalB, long& ccLag) {

	double aAvg = 0, bAvg = 0;
	double aSqrSum = 0, bSqrSum = 0;
	double num = 0, denom = 1;
	int sampleNr = signalA.size();

	for (int i = 0; i < sampleNr; i++) {
		aAvg += signalA[i];
		bAvg += signalB[i];
	}
	aAvg /= sampleNr;
	bAvg /= sampleNr;


	for (int i = 0; i < sampleNr; i++) {
		if (i + ccLag < 0)
			continue;
		if (i + ccLag >= sampleNr)
			break;
		num += (signalA[i] - aAvg)*(signalB[i + ccLag] - bAvg);
		aSqrSum += pow(signalA[i] - aAvg, 2);
		bSqrSum += pow(signalB[i + ccLag] - bAvg, 2);
	}
	denom = sqrt(aSqrSum*bSqrSum);
	if (denom == 0) {
		return 0;
	}
	return num / denom;
}

vectorStats SignalProcessing::vectorStatistics(std::vector<double> vec) {

	vectorStats stats;
	int vecSize = vec.size();
	double maxVal = *std::max_element(vec.begin(), vec.end());

	if (vecSize == 0) return stats;

	std::sort(vec.begin(), vec.end());
	stats.min = vec[0];
	stats.max = vec[vecSize - 1];

	if (vec.size() < 5) {
		for (unsigned i = 0; i < vec.size(); ++i) {
			stats.median += vec[vec.size() / 2];
		}
		stats.median /= vec.size();
		stats.firstQ = stats.min;
		stats.thirdQ = stats.max;
		return stats;
	}


	std::vector<double> firstQ_Vec;
	std::vector<double> thirdQ_Vec;

	double medianIdxInt;
	double medianIdx = vectorMedian(vec, stats.median);

	if (std::modf(medianIdx, &medianIdxInt) == 0.0) {	// median is a datum
		for (int i = 0; i <= medianIdxInt; ++i)
			firstQ_Vec.push_back(vec[i]);
		for (int i = medianIdxInt; i < vecSize; ++i)
			thirdQ_Vec.push_back(vec[i]);
		vectorMedian(firstQ_Vec, stats.firstQ);
		vectorMedian(thirdQ_Vec, stats.thirdQ);
	}
	else {	// median is not a datum
		for (int i = 0; i <= medianIdxInt + 1; ++i)
			firstQ_Vec.push_back(vec[i]);
		for (int i = medianIdxInt + 1; i < vecSize + 1; ++i)
			thirdQ_Vec.push_back(vec[i]);
		vectorMedian(firstQ_Vec, stats.firstQ);
		vectorMedian(thirdQ_Vec, stats.thirdQ);
	}

	// the ends of the whiskers represent the lowest datum still within 1.5 IQR of the lower quartile, and the highest datum still within 1.5 IQR of the upper quartile (Tukey boxplot)
	double iqr = stats.thirdQ - stats.firstQ;
	for (int i = vec.size() - 1; i > 0; --i) {
		if (vec[i] <= stats.thirdQ + 1.5*iqr) {
			stats.max = vec[i];
			break;
		}
	}
	for (int i = 0; i < vec.size(); ++i) {
		if (vec[i] >= stats.firstQ - 1.5*iqr) {
			stats.min = vec[i];
			break;
		}
	}

	if (vecSize < 2)
		stats.variance = 0;
	else {
		for (int i = 0; i < vecSize; ++i) {
			stats.mean += vec[i];
		}
		stats.mean /= vecSize;

		for (int i = 0; i < vecSize; ++i) {
			stats.variance += ((vec[i] - stats.mean)*(vec[i] - stats.mean));
		}
		stats.variance /= (vecSize - 1);
	}

	stats.stdrdDev = sqrt(stats.variance);

	return stats;

}

double SignalProcessing::vectorMedian(std::vector<double> vec, double &median) {

	double medianIdx = 0;

	if (vec.empty()) return 0;

	/*std::string vect;
	for (int i = 0; i < vec.size(); ++i) {
	vect += std::to_string(vec[i]) + ", ";
	}*/

	int vecSize = vec.size();
	if (vec.size() % 2 == 0) {
		medianIdx = ((vecSize / 2 - 1) + (vecSize / 2)) / 2;
		median = (vec[vecSize / 2 - 1] + vec[vecSize / 2]) / 2;
	}
	else {
		medianIdx = vecSize / 2;
		median = vec[vecSize / 2];
	}

	return medianIdx;
}

bool SignalProcessing::characterizeBandpassRippleFilter(double fOne, double fTwo) {

	_mkdir("FilterCharacterization/");
	
	unsigned firOrder = FIR_RIPPLE_ORDER;
	double m_samplingRate = 2048;
	unsigned durationSamples = 2*m_samplingRate;	//2 seconds

	std::vector<double> frequencyVec;
	std::vector<double> gainVec, delayDegrees, delayRadians;
	double timePrec = 1.0 / m_samplingRate;
	double freqPrec = 1;
	double f_1 = 0, f_0 = 0, f_2 = 0;
	double pi = 3.14159265359;
	//Transfer Function
	for (double freqDegrees = 1; freqDegrees < 800; freqDegrees += freqPrec) {
		int sineFirstPeak_sample = 0;
		int filteredFirstPeak_sample = 0;
		double max = 0;

		matrixStd sineMatrix(1);
		sineMatrix[0].resize(durationSamples);
		std::vector<double> sineVec;
													// generate the sine signal for the given freqDegrees
		for (int i = 0; i < durationSamples; ++i) {
			double sampleTime = ((double)i)*timePrec;
			double sinVal = sin(2.0 * pi *freqDegrees * sampleTime);
			sineMatrix[0][i] = (sinVal);
			sineVec.push_back(sinVal);

			if (sineVec[i]> max) {
				max = sineVec[i];
				sineFirstPeak_sample = i;
			}

		}


		//Filter in ripple range	
		matrixStd filtOutMatrix, notchedHighPassed;
		std::vector<double> filtOutVec;
		rippleBP(sineMatrix, filtOutMatrix, fOne, fTwo, m_samplingRate);	//Ripple_Range
		max = 0;
		for (int i = 0; i < filtOutMatrix[0].size(); ++i) {
			double filtVal = filtOutMatrix[0][i];
			filtOutVec.push_back(filtVal);


			if (filtOutVec[i]> max) {
				max = filtOutVec[i];
				filteredFirstPeak_sample = i;
			}
		}

		double maxAmplitudeSine = *std::max_element(sineVec.begin(), sineVec.end()) - *std::min_element(sineVec.begin(), sineVec.end());
		double maxAmplitudeFiltOut = *std::max_element(filtOutVec.begin(), filtOutVec.end()) - *std::min_element(filtOutVec.begin(), filtOutVec.end());

		double gain_dB = maxAmplitudeSine > 0 ? 20 * log10(maxAmplitudeFiltOut / maxAmplitudeSine) : 0;
		//gain_dB = (maxAmplitudeFiltOut / maxAmplitudeSine);

		double delayTime = abs(((double)sineFirstPeak_sample / m_samplingRate) - ((double)filteredFirstPeak_sample / m_samplingRate)) * -1;
		double phaseDelayDegrees = 360 * freqDegrees * delayTime;
		double phaseDelayRadians = 2 * pi * freqDegrees * delayTime;

		gainVec.push_back(gain_dB);
		delayDegrees.push_back(phaseDelayDegrees);
		delayRadians.push_back(phaseDelayRadians);
		frequencyVec.push_back(freqDegrees);
	}
	std::string signalName = "FilterCharacterization/Band-pass, FreqResponse," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::plotXY(frequencyVec, gainVec, signalName, 0, 0);
	Plotter::writeVectorsToFile(frequencyVec, gainVec, signalName);

	signalName = "FilterCharacterization/Band-pass, PhaseShift, Degrees," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayDegrees, signalName);

	signalName = "FilterCharacterization/Band-pass, PhaseShift, Radians," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayRadians, signalName);

	int impulseResponseLength = 2*m_samplingRate;
	std::vector<double> time;
	matrixStd impulse(1), impulseResponse;
	for (int i = 0; i < impulseResponseLength; ++i) {
		time.push_back(i / m_samplingRate);
		impulse[0].push_back(0);
	}
	impulse[0][impulseResponseLength / 2] = 1;

	rippleBP(impulse, impulseResponse, fOne, fTwo, m_samplingRate);	//Ripple_Range

	double peaksRateNr = 0;
	for (int i = 1; i < impulseResponseLength - 1; ++i) {
		double val = impulseResponse[0][i];

		double before = impulseResponse[0][i - 1];
		double center = val;
		double after = impulseResponse[0][i + 1];
		if (compareDoublesLarger(center, before) && compareDoublesLarger(center, after) || compareDoublesSmaller(center, before) && compareDoublesSmaller(center, after)) {
			peaksRateNr++;
		}
	}

	std::string impusleResponseName = "FilterCharacterization/Band-pass, ImpulseResponse" + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder) + ", Oscillations" + std::to_string(peaksRateNr/2);
	Plotter::plotXY(time, impulseResponse[0], impusleResponseName, 0, 0);
	Plotter::writeVectorsToFile(time, impulseResponse[0], impusleResponseName);
	
	return true;
}

bool SignalProcessing::characterizeBandpassFastRippleFilter(double fOne, double fTwo) {

	_mkdir("FilterCharacterization/");

	unsigned firOrder = FIR_FR_ORDER;
	double m_samplingRate = 2048;
	unsigned durationSamples = 2 * m_samplingRate;	//2 seconds

	std::vector<double> frequencyVec;
	std::vector<double> gainVec, delayDegrees, delayRadians;
	double timePrec = 1.0 / m_samplingRate;
	double freqPrec = 1;
	double f_1 = 0, f_0 = 0, f_2 = 0;
	double pi = 3.14159265359;
	//Transfer Function
	for (double freqDegrees = 1; freqDegrees < 800; freqDegrees += freqPrec) {
		int sineFirstZC_sample = 0;
		int filteredFirstZC_sample = 0;
		matrixStd sineMatrix(1);
		sineMatrix[0].resize(durationSamples);
		std::vector<double> sineVec;
		// generate the sine signal for the given freqDegrees
		for (int i = 0; i < durationSamples; ++i) {
			double sampleTime = ((double)i)*timePrec;
			double sinVal = sin(2.0 * pi *freqDegrees * sampleTime);
			sineMatrix[0][i] = (sinVal);
			sineVec.push_back(sinVal);
			if (i > 0 && sineFirstZC_sample == 0) {
				bool posZC = sineVec[i - 1] >= 0 && sineVec[i] < 0;
				bool negZC = sineVec[i - 1] <= 0 && sineVec[i] > 0;
				if (posZC || negZC)
					sineFirstZC_sample = i;
			}
		}


		//Filter in ripple range	
		matrixStd filtOutMatrix, notchedHighPassed;
		std::vector<double> filtOutVec;
		fastRippleBP(sineMatrix, filtOutMatrix, fOne, fTwo, m_samplingRate);	//Ripple_Range
		for (int i = 0; i < filtOutMatrix[0].size(); ++i) {
			double filtVal = filtOutMatrix[0][i];
			filtOutVec.push_back(filtVal);

			if (i > 0 && filteredFirstZC_sample == 0) {
				bool posZC = filtOutVec[i - 1] >= 0 && filtOutVec[i] < 0;
				bool negZC = filtOutVec[i - 1] <= 0 && filtOutVec[i] > 0;
				if (posZC || negZC)
					filteredFirstZC_sample = i;
			}
		}

		double maxAmplitudeSine = *std::max_element(sineVec.begin(), sineVec.end()) - *std::min_element(sineVec.begin(), sineVec.end());
		double maxAmplitudeFiltOut = *std::max_element(filtOutVec.begin(), filtOutVec.end()) - *std::min_element(filtOutVec.begin(), filtOutVec.end());

		double gain_dB = maxAmplitudeSine > 0 ? 20 * log10(maxAmplitudeFiltOut / maxAmplitudeSine) : 0;
		//gain_dB = (maxAmplitudeFiltOut / maxAmplitudeSine);

		double delayTime = ((double)filteredFirstZC_sample / m_samplingRate) - ((double)sineFirstZC_sample / m_samplingRate);
		double phaseDelayDegrees = 360 * freqDegrees * delayTime;
		double phaseDelayRadians = 2 * pi * freqDegrees * delayTime;

		gainVec.push_back(gain_dB);
		delayDegrees.push_back(phaseDelayDegrees);
		delayRadians.push_back(phaseDelayRadians);
		frequencyVec.push_back(freqDegrees);
	}
	std::string signalName = "FilterCharacterization/Band-pass, FreqResponse," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::plotXY(frequencyVec, gainVec, signalName, 0, 0);
	Plotter::writeVectorsToFile(frequencyVec, gainVec, signalName);

	signalName = "FilterCharacterization/Band-pass, PhaseShift, Degrees," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayDegrees, signalName);

	signalName = "FilterCharacterization/Band-pass, PhaseShift, Radians," + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayRadians, signalName);

	int impulseResponseLength = 2 * m_samplingRate;
	std::vector<double> time;
	matrixStd impulse(1), impulseResponse;
	for (int i = 0; i < impulseResponseLength; ++i) {
		time.push_back(i / m_samplingRate);
		impulse[0].push_back(0);
	}
	impulse[0][impulseResponseLength / 2] = 1;

	fastRippleBP(impulse, impulseResponse, fOne, fTwo, m_samplingRate);	//Ripple_Range

	double peaksRateNr = 0;
	for (int i = 1; i < impulseResponseLength - 1; ++i) {
		double val = impulseResponse[0][i];

		double before = impulseResponse[0][i - 1];
		double center = val;
		double after = impulseResponse[0][i + 1];
		if (compareDoublesLarger(center, before) && compareDoublesLarger(center, after) || compareDoublesSmaller(center, before) && compareDoublesSmaller(center, after)) {
			peaksRateNr++;
		}
	}

	std::string impusleResponseName = "FilterCharacterization/Band-pass, ImpulseResponse" + std::to_string((int)fOne) + "-" + std::to_string((int)fTwo) + ", FIR Order-" + std::to_string(firOrder) + ", Oscillations" + std::to_string(peaksRateNr / 2);
	Plotter::plotXY(time, impulseResponse[0], impusleResponseName, 0, 0);
	Plotter::writeVectorsToFile(time, impulseResponse[0], impusleResponseName);

	return true;
}

bool SignalProcessing::characterizeHighpassFilter(double cuttofFreq) {

	_mkdir("FilterCharacterization/");

	unsigned firOrder = FIR_HP_ORDER;
	double m_samplingRate = 2048;
	unsigned durationSamples = 2 * m_samplingRate;	//2 seconds

	std::vector<double> frequencyVec;
	std::vector<double> gainVec, delayDegrees, delayRadians;
	double timePrec = 1.0 / m_samplingRate;
	double freqPrec = 1;
	double f_1 = 0, f_0 = 0, f_2 = 0;
	double pi = 3.14159265359;
	//Transfer Function
	for (double freqDegrees = 1; freqDegrees < 800; freqDegrees += freqPrec) {
		int sineFirstZC_sample = 0;
		int filteredFirstZC_sample = 0;

		matrixStd sineMatrix(1);
		sineMatrix[0].resize(durationSamples);
		std::vector<double> sineVec;
		// generate the sine signal for the given freqDegrees
		for (int i = 0; i < durationSamples; ++i) {
			double sampleTime = ((double)i)*timePrec;
			double sinVal = sin(2.0 * pi *freqDegrees * sampleTime);
			sineMatrix[0][i] = (sinVal);
			sineVec.push_back(sinVal);
			if (i > 0 && sineFirstZC_sample == 0) {
				bool posZC = sineVec[i - 1] >= 0 && sineVec[i] < 0;
				bool negZC = sineVec[i - 1] <= 0 && sineVec[i] > 0;
				if(posZC || negZC)
					sineFirstZC_sample = i;
			}
		}


		//Filter in ripple range	
		matrixStd filtOutMatrix, notchedHighPassed;
		std::vector<double> filtOutVec;
		applyHighPassFilterToSignal(sineMatrix, filtOutMatrix, cuttofFreq, m_samplingRate);		//Highpass_1Hz	
		for (int i = 0; i < filtOutMatrix[0].size(); ++i) {
			double filtVal = filtOutMatrix[0][i];
			filtOutVec.push_back(filtVal);

			if (i > 0 && filteredFirstZC_sample == 0) {
				bool posZC = filtOutVec[i - 1] >= 0 && filtOutVec[i] < 0;
				bool negZC = filtOutVec[i - 1] <= 0 && filtOutVec[i] > 0;
				if (posZC || negZC)
					filteredFirstZC_sample = i;
			}
		}

		double maxAmplitudeSine = *std::max_element(sineVec.begin(), sineVec.end()) - *std::min_element(sineVec.begin(), sineVec.end());
		double maxAmplitudeFiltOut = *std::max_element(filtOutVec.begin(), filtOutVec.end()) - *std::min_element(filtOutVec.begin(), filtOutVec.end());

		double gain_dB = maxAmplitudeSine > 0 ? 20 * log10(maxAmplitudeFiltOut / maxAmplitudeSine) : 0;
		//gain_dB = (maxAmplitudeFiltOut / maxAmplitudeSine);

		double delayTime = ((double)filteredFirstZC_sample / m_samplingRate) - ((double)sineFirstZC_sample / m_samplingRate);
		double phaseDelayDegrees = 360 * freqDegrees * delayTime;
		double phaseDelayRadians = 2 * pi * freqDegrees * delayTime;

		gainVec.push_back(gain_dB);
		delayDegrees.push_back(phaseDelayDegrees);
		delayRadians.push_back(phaseDelayRadians);
		frequencyVec.push_back(freqDegrees);
	}
	std::string signalName = "FilterCharacterization/High-pass, FreqResponse," + std::to_string((int)cuttofFreq) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::plotXY(frequencyVec, gainVec, signalName, 0, 0);
	Plotter::writeVectorsToFile(frequencyVec, gainVec, signalName);

	signalName = "FilterCharacterization/High-pass, PhaseShift, Degrees," + std::to_string((int)cuttofFreq) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayDegrees, signalName);
	
	signalName = "FilterCharacterization/High-pass, PhaseShift, Radians," + std::to_string((int)cuttofFreq) + ", FIR Order-" + std::to_string(firOrder);
	Plotter::writeVectorsToFile(frequencyVec, delayRadians, signalName);

	int impulseResponseLength = 2 * m_samplingRate;
	std::vector<double> time;
	matrixStd impulse(1), impulseResponse;
	for (int i = 0; i < impulseResponseLength; ++i) {
		time.push_back(i / m_samplingRate);
		impulse[0].push_back(0);
	}
	impulse[0][impulseResponseLength / 2] = 1;

	applyHighPassFilterToSignal(impulse, impulseResponse, cuttofFreq, m_samplingRate);	//Ripple_Range

	double peaksRateNr = 0;
	for (int i = 1; i < impulseResponseLength - 1; ++i) {
		double val = impulseResponse[0][i];

		double before = impulseResponse[0][i - 1];
		double center = val;
		double after = impulseResponse[0][i + 1];
		if (compareDoublesLarger(center, before) && compareDoublesLarger(center, after) || compareDoublesSmaller(center, before) && compareDoublesSmaller(center, after)) {
			peaksRateNr++;
		}
	}

	std::string impusleResponseName = "FilterCharacterization/High-pass, ImpulseResponse" + std::to_string((int)cuttofFreq) + ", FIR Order-" + std::to_string(firOrder) + ", Oscillations" + std::to_string(peaksRateNr / 2);
	Plotter::plotXY(time, impulseResponse[0], impusleResponseName, 0, 0);
	Plotter::writeVectorsToFile(time, impulseResponse[0], impusleResponseName);

	return true;
}