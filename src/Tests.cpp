#include <random>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <direct.h>

#include "SignalProcessing.h"
#include "Tests.h"

#define RPPL_LOW_FREQ			80
#define RPPL_HIGH_FREQ			250
#define FAST_RPPL_LOW_FREQ		250
#define FAST_RPPL_HIGH_FREQ		500	

void characterizeFilter(std::string outputDirectory, double samplingRate, double firOrder) {

	double m_samplingRate = 2048;// samplingRate;
	std::string m_outputDirectory = outputDirectory;
	double m_firOrder = firOrder;

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	mkdir(outputPath.c_str());

	int signalLength = 0.5 * m_samplingRate;
	matrixStd impulse(1);
	matrixStd impulseResponseRipple, impulseResponseFR;

	for (unsigned i = 0; i < signalLength; ++i) {
		impulse[0].push_back(0);
	}
	impulse[0][signalLength / 2] = 10;

	// Filter in the frequency band of EOI
	SignalProcessing::applyBandPassFilterToSignal(impulse, impulseResponseRipple, RPPL_LOW_FREQ, RPPL_HIGH_FREQ, m_samplingRate);
	SignalProcessing::applyBandPassFilterToSignal(impulse, impulseResponseFR, FAST_RPPL_LOW_FREQ, FAST_RPPL_HIGH_FREQ, m_samplingRate);

	std::ofstream impRespFile;
	std::string filename = m_outputDirectory + "/MOSSDET_Output/ImpulseResponse,FIR Order-" + std::to_string(m_firOrder) + ".xls";
	remove(filename.c_str());
	impRespFile.open(filename);
	impRespFile << "Time\tImpulse\tImpulseResponseRipple\tImpulseResponseFastRipple\n";
	for (int i = 0; i < impulseResponseRipple[0].size(); ++i) {
		impRespFile << i / m_samplingRate  << "\t" << impulse[0][i] << "\t" << impulseResponseRipple[0][i] << "\t" << impulseResponseFR[0][i] << "\n";
	}
	impRespFile.close();
}