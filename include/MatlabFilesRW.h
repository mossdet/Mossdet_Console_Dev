#pragma once

#include "DataTypes.h"


class MatlabFilesRW {

public:
	MatlabFilesRW(std::string filename, double samplingRate, std::string outputPath, bool verbose);
	int readMatlabSamples(long startSample, long nrSamplesToRead, std::vector<std::vector<double>> &m_reformattedSignalSamples);
	int writeMatlabSamples(std::vector<EOI_Event> &allDetections);
	long getSampleCount();

private:
	std::string m_matFilename;
	std::string m_samplingRate;
	std::string m_outputPath;
	bool m_verbose;
};
