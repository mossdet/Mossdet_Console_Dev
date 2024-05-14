#pragma once

#include "DataTypes.h"

/*#define FIR_RIPPLE_ORDER 18
#define FIR_FR_ORDER 10
#define FIR_HP_ORDER 64*/

#define FIR_RIPPLE_ORDER 256
#define FIR_FR_ORDER 256
#define FIR_HP_ORDER 256

/*#define FIR_RIPPLE_ORDER 64
#define FIR_FR_ORDER 64
#define FIR_HP_ORDER 64*/

#define FIR_ORDER FIR_HP_ORDER//SS-512; HFO-256
class SignalProcessing
{
public:
	SignalProcessing();
	~SignalProcessing();

	static bool applyBandPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate);
	static bool applyHighPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesHP, double dLowCutOffFreq, double samplingRate);
	static bool applyLowPassFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesHP, double dLowCutOffFreq, double samplingRate);
	static bool applyBandStopFilterToSignal(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBS, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate);

	static bool rippleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate);
	static bool fastRippleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate);
	static bool spindleBP(matrixStd& signalSamples, matrixStd& filteredSignalSamplesBP, double dLowCutOffFreq, double dHighCutOffFreq, double samplingRate);

	static bool firstDerivative(std::vector<double>& signal, std::vector<double>& firstDerivative, double startSample, double sampleLength, double samplingRate);

	static unsigned morletWaveletTransform(matrixStd &inData, double dSamplingFreq, double dStartFreq, double dEndFreq, double dDeltaFreq, int waveletOscills, matrixStd& waveletOutput);
	static void getMatrixPower(matrixStd& matrixOutput);

	bool setVarianceData(const std::vector<double> data);
	double runVariance(long startSample, long sampleLength);
	static double getPearsonCorrCoeff(std::vector<double>& signalA, std::vector<double>& signalB, long& ccLag);

	static vectorStats vectorStatistics(std::vector<double> vec);
	static double vectorMedian(std::vector<double> vec, double &median);
	static bool characterizeBandpassRippleFilter(double fOne, double fTwo);
	static bool characterizeBandpassFastRippleFilter(double fOne, double fTwo);
	static bool characterizeHighpassFilter(double cuttofFreq);

	static double getOrder() {
		return FIR_ORDER;
	}

	static bool compareDoublesEquals(double dFirstVal, double dSecondVal) {
		return std::fabs(dFirstVal - dSecondVal) < std::numeric_limits<double>::epsilon() * 100;
	}

	static bool compareDoublesLarger(double dFirstVal, double dSecondVal) {
		return dFirstVal - dSecondVal > std::numeric_limits<double>::epsilon() * 100;
	}

	static bool compareDoublesSmaller(double dSmaller, double dBigger) {
		return dBigger - dSmaller > std::numeric_limits<double>::epsilon() * 100;
	}


protected:
	std::vector<double> m_varData;

	//static const double FIR_ORDER = 256;//SS-512;

	unsigned int bandpassRippleOrder;
	unsigned int bandpassFastRippleOrder;
	unsigned int highPassOrder;

};

