#pragma once

#include <string>
#include <math.h>

#include "DataTypes.h"

#define REJECT_ARTIFACTS
#ifdef REJECT_ARTIFACTS
	#undef REJECT_ARTIFACTS
#endif

class SVM_Detector {
public:
	SVM_Detector(std::string eoiName, std::string patientName, std::string functionName, std::string trainingDataFileName, 
		std::vector<unsigned> selectedFeatures, std::string outputPath, bool saveTextFile);
	~SVM_Detector();

	bool detectEpochs(double chNr);
	bool getEEG_EventsFromDetectedEpochs(unsigned sharedChs, std::string channName);
	bool generateDetectionsAndFeaturesFile(double onsetMinDiff, double minSharedPercent);

	bool compareDoublesEquals(double dFirstVal, double dSecondVal) {
		return std::fabs(dFirstVal - dSecondVal) < std::numeric_limits<double>::epsilon() * 100;
	}

	bool compareDoublesLarger(double dFirstVal, double dSecondVal) {
		return dFirstVal - dSecondVal > std::numeric_limits<double>::epsilon() * 100;
	}

	bool compareDoublesSmaller(double dSmaller, double dBigger) {
		return dBigger - dSmaller > std::numeric_limits<double>::epsilon() * 100;
	}

	void setInputEvents(std::vector<EOI_Epoch>  &inputEpochs) {
		m_inputEvents = inputEpochs;
		m_signalDuration = m_inputEvents.back().time - m_inputEvents[0].time;
	}

	void setSamplingRate(double samplingRate) {
		m_samplingRate = samplingRate;
	}

	void setGammaAndNu(double gamma, double nu) {
		m_gamma = gamma;
		m_nu = nu;
	}

	void clearReadVals() {
		m_inputEvents.clear();
	}

	void initEpochVectors(unsigned nrChannels) {
		m_chSpecifificDetecEpochs.resize(nrChannels);
	}

	void setProbThreshold(double varMedian, double val) {
		m_probThreshold = val;
	}

	void clearDetectedEpochs() {
		m_chSpecifificDetecEpochs.clear();
	}

	void clearEpochVectors() {
		m_chSpecifificDetecEpochs.clear();
		m_detectedEOI.clear();
	}

	void setFunctionName(std::string functionName) {
		m_savedFunctionName = functionName;
	}

	void setTrainingFileName(std::string trainingFileName) {
		m_trainingDataFileName.clear();
		m_trainingDataFileName = trainingFileName;
	}

	double getAvg(std::vector<double> &featVec) {
		double avg = 0;
		for (size_t i = 0; i < featVec.size(); ++i)
			avg += featVec[i];
		return avg / (double)featVec.size();
	}

	bool fileExists(const char *fileName);

	void setEOI_Name(std::string eoiName) {
		m_eoiName = eoiName;
	}

	void setCI(bool ciEnable) {
		m_useConfidenceInterval = ciEnable;
	}

	void setPatientPath(std::string patientPath) {
		m_patientPath = patientPath;
	}

	void setSelectedFeatures(std::vector<unsigned> selectedFeats) {
		m_selectedFeats = selectedFeats;
	}

	void setAbsoluteSignalLength(double absSigLength) {
		m_absoluteSignalLength = absSigLength;
	}

	bool isChannelAnnotated() {

		//bool channelisAnnotated = m_visualEOI_Joined.size() != 0;
		bool channelisAnnotated = m_visualEOI_Joined.size() > 0;

		return  channelisAnnotated;
	}

	bool getEvents(std::vector<EOI_Event> &readEvents);

	void setFileStartDateAndTime (DateAndTime fileStartDateAndTime){ m_fileStartDateAndTime  = fileStartDateAndTime; }

private:
	// Classifier parameters
	std::string m_trainingDataFileName;
	std::vector<unsigned> m_selectedFeats;
	double m_gamma;
	double m_nu;
	double m_probThreshold;
	double m_epochLength;
	std::string m_savedFunctionName;
	std::string m_eoiName;

	// Characterized Signal parameters
	std::string m_patientName;
	std::string m_patientPath;
	std::string m_filename;
	double m_avgMarkDuration;
	double m_expertAgreement;
	double m_interExpertKappa;
	std::vector<EOI_Epoch> m_inputEvents;
	unsigned m_nrEvents;
	double m_samplingRate;
	double m_signalDuration;
	double m_epochsStartTime;
	double m_epochsEndTime;
	double m_absoluteSignalLength;
	unsigned m_nrSpikes;
	unsigned m_nrSpikes_wo_HFO;

	// Detections Results
	unsigned m_nrShort_EOI;
	unsigned m_nrLong_EOI;
	unsigned m_nrNoOscills_EOI;
	unsigned m_nrNoHistoryCompliant_EOI;
	unsigned m_artefactualDetections;
	unsigned m_nrSpikesCausingFalseHFO;
	double m_avgDetectionDuration;
	std::vector<std::vector<bool>> m_chSpecifificDetecEpochs;
	std::vector<EOI_Event> m_detectedEOI;
	std::vector<EOI_Event> m_visualEOI;
	std::vector<EOI_Event> m_visualEOI_Joined;
	
	std::vector<EOI_Event> m_visualEventsExpOne;
	std::vector<EOI_Event> m_visualEventsExpTwo;
	std::vector<EOI_Event> m_visualSpikes;
	std::vector<EOI_Event> m_visualSpikesWithHFO;
	std::vector<EOI_Event> m_visualSpikes_wo_HFO;

	std::string m_expOneName;
	std::string m_expTwoName;


	//Initialize File Paths
	std::string m_outputPath;
	std::string m_patPath;
	std::string m_patDetectionFilesPath;

	Performance m_svmPerformance;

	bool m_useConfidenceInterval;

	std::vector<bool> m_spikeDetectedEpochs;

	DateAndTime m_fileStartDateAndTime;
	bool m_saveTextFile;
};