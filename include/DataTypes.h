#pragma once

#include <vector>

#define _T(x)               x 

typedef std::vector<std::vector<double>> matrixStd;	//(rows, columns)

struct Channel {
	unsigned channelNr = 0;
	std::vector<double> data;
};

struct ContactNames {
	unsigned contactGlobalIdx = 0;
	unsigned contactNr = 0;
	std::string electrodeName;
	std::string contactName;

};

struct MontageNames {
	std::string montageName;
	unsigned firstContactGlobalIdx = 0, secondContactGlobalIdx;
	unsigned firstContactNr = 0, secondContactNr;
	int montageMOSSDET_Nr;
	std::string firstElectrodeName, secondElectrodeName;
};

struct EOI_Event {
	std::string channelName;
	std::string lowNeighborChannelName;
	std::string highNeighborChannelName;
	std::string description;
	unsigned channelNr = 0;
	unsigned lowNeighborChannelNr = 0;
	unsigned highNeighborChannelNr = 0;
	double startTime = 0;
	double endTime = 0;
	double duration = 0;
	double label = 0;
	double amplitude = 0;
	double power = 0;
	double spectralPeak = 0;
	double backgroundAmplitude = 0;
	double backgroundPower = 0;
	double backgroundStdDev = 0;
};

struct EOI_Epoch {
	unsigned channel;
	double time;
	std::vector<double> inputVals;
	std::vector<double> targetVals;
	std::vector<double> classifierOutput;
};

struct Performance {
	double truePositives = 0, falsePositives = 0, trueNegatives = 0, falseNegatives = 0;
};

struct DateAndTime{
	long day = 0, month = 0, year = 0;
	double hours = 0, minutes= 0, seconds= 0, milliSeconds = 0;
};

class Patient {

public:
	Patient() {
		startHour = 0.0; startMin = 0.0; startSec = 0.0;
	}

	std::string patientName;
	std::vector<std::string> DataFilenames;
	std::string TrainingDataFilename;
	std::string MarkingsFn;
	std::string RefinedMarkingsFn;
	std::string selectedFeaturesFnSCH;
	std::string selectedFeaturesFnMCH;

	std::vector<std::pair <int, std::string>> montageChanLabels;

	unsigned startHour, startMin;
	double startSec;
	double samplingRate;

	double totalStartSecs() {
		double totStartSecs = startHour * 3600 + startMin * 60 + startSec;
		return totStartSecs;
	}

	void clear() {
		patientName.clear();
		DataFilenames.clear();
		MarkingsFn.clear();
		RefinedMarkingsFn.clear();
		startHour = startMin = startSec = 0;
	}

};

struct vectorStats {
	double median = 0;
	double firstQ = 0;
	double thirdQ = 0;
	double max = 0;
	double min = 0;
	double mean = 0;
	double variance = 0;
	double stdrdDev = 0;
};
