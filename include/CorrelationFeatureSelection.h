#pragma once

#include <string>
#include <math.h>

#include "HFO_EventsHandler.h"
#include "DataTypes.h"

class CorrelationFeatureSelection {
public:
	CorrelationFeatureSelection(unsigned nrChannels, std::string patientPath, std::string patientName);

	bool getEOI_MarkFeatureCorrelation(bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes, EOI_Features &epochFeatures);
	bool getCFS(std::vector<double> featMarkCorrChannAvg, std::vector<std::vector<double>> featFeatCorrChannAvg, bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes);
	long getMaxCrossCorrelationPearson(std::vector<double> &A, std::vector<double> &B, double& corrCoeff, long maxNegLag, long maxPosLag);
	long getMaxCrossCorrelationRho(std::vector<double> &signalA, std::vector<double> &signalB, long startSample, long length, double& corrCoeff, long maxNegLag, long maxPosLag);
	
	bool relabelDataHFO(bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, EOI_Features &epochFeatures);
	bool relabelDataSpikesHFO(EOI_Features &epochFeatures);
	bool relabelDataSpikesAlone(EOI_Features &epochFeatures);
	bool relabelDataAllSpikes(EOI_Features &epochFeatures);

	bool vectorizeChannels(EOI_Features &epochFeatures, EOI_Features &vectorizedFeatureChannels);

	void setMultichannel(bool multichannel) {
		m_multichannel = multichannel;
	}

	bool compareDoubles(double dFirstVal, double dSecondVal){
		return (dFirstVal < dSecondVal + 0.1) && (dFirstVal > dSecondVal - 0.1);
		//return std::fabs(dFirstVal - dSecondVal) < std::numeric_limits<double>::epsilon()* 100;
	}

	void setFeatsToExclude(std::vector<int> featsToExclude){
		m_featsToExclude = featsToExclude;
	}

private:

	unsigned m_nrChannels;
	bool m_multichannel;
	std::string m_patientName;
	std::string m_patientPath;

	std::vector<int> m_featsToExclude;
};