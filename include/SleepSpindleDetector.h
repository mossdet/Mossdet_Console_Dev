#pragma once

#include <string>
#include "DataTypes.h"
#include "SS_EventsHandler.h"

typedef struct WaveletParametersSS {
	long lOverlapSamples = 0;
	unsigned int waveletOscillations = 6;
	double dStartFrequency = 0;
	double dEndFrequency = 0;
	double dDeltaFrequency = 0;
	unsigned nrComponents;
};

class SleepSpindleDetector
{
public:

	SleepSpindleDetector(std::string strInputFileName, std::string detectorPath, std::string strDirName, int patIdx, int montage, std::vector<std::string>* selectedUnipolarChanels, std::vector<std::string>* selectedBipolarChanels,
		std::string& eoiName, bool useChannsFileBool, double startTime, double endTime, bool compareDetectionsWithVisalMarks, double samplingRate, bool verbose, bool saveTextFile);

	long readEEG();
	bool characterizeEEG();
	bool detect(SS_Features& signalFeatures, long firstLocalSampleToRead);
	bool generateEpochCharacterizationSpindles(SS_Features& epochFeatures);


	int getWaveletCausedDelay(void);
	int getFilteringCausedDelay(void);
	bool getFeatureAnalysisAllSignals(matrixStd& rawSignal, SS_Features& epochFeatures, long rewind, long overlap, long firstSampleToRead, long& samplesToRead, long epochLength, int signalType);
	void getSignalParams(double& samplingRate, double& analysisDuration, int& nrChannels);

	bool getWaveletFrequencyeAnalysisAllSignal(matrixStd& rawSignal, SS_Features& epochFeatures, long rewind, long overlap, long firstSampleToRead, long& samplesToRead,
		int firstChannel, int lastChannel, long epochLength);

	std::string getPatientNameAndPath();

	bool saveSelectedEEG_DataToFile(long firstSampleToRead, unsigned rewind, unsigned overlap, matrixStd& signal, matrixStd& notchedSignal, matrixStd& bpSignal, unsigned channelToSave);
	bool plotSelectedEE_Channel(long firstSampleToRead, unsigned rewind, unsigned overlap, matrixStd& signal, matrixStd& notchedSignal, matrixStd& rippleSignal, matrixStd& fastRippleSignal, unsigned channelToSave);
	void writeVectorToFile(int channNr, double time, std::vector<double> frequency, std::vector<double> vectorB);

	void readChannelsFromFile(void);
	void getSpecificStartAndEndTimes();
	void readAndWriteEOI_VisualAnnotations();
	bool getChannelNumber(std::string markChannelName, unsigned int& channNr);
	bool getNeighboringNumber(long annotChannNr, long& lowerNeighborChNr, long& upperNeighborChNr);
	void generate_EOI_ChannelsFile(void);

	std::string ExePath();

	bool getPerformance() {
		return m_compareDetectionsWithVisalMarks;
	}

	bool compareDoublesEquals(double dFirstVal, double dSecondVal) {
		return (dFirstVal < dSecondVal + 0.001) && (dFirstVal > dSecondVal - 0.001);
		//return std::fabs(dFirstVal - dSecondVal) < std::numeric_limits<double>::epsilon()* 100;
	}

	bool compareDoublesLarger(double dFirstVal, double dSecondVal) {
		return dFirstVal - dSecondVal > std::numeric_limits<double>::epsilon() * 100;
	}

	bool compareDoublesSmaller(double dSmaller, double dBigger) {
		return dBigger - dSmaller > std::numeric_limits<double>::epsilon() * 100;
	}

	bool fileExists(const char* fileName) {
		std::ifstream infile(fileName);
		return infile.good();
	}

	//private:
	matrixStd m_reformattedSignalSamples;

	//private:
	std::string m_strFileName;
	std::string m_patientName;
	std::string m_patientPath;
	std::string m_exePath;
	std::string m_outputDirectory;
	std::string m_marksFolder;
	std::string m_eoiName;

	std::vector<ContactNames> m_unipolarLabels;
	std::vector<MontageNames> m_montageLabels;

	std::vector<std::string> m_selectedUnipolarChanels;
	std::vector<std::string> m_selectedBipolarChanels;

	std::vector<std::vector<EOI_Event>> m_markedSS;

	long m_totalSampleCount;
	double m_samplingRate;
	int m_nrChannels;
	long long m_fileAnalysisStartSample;
	long long m_fileAnalysisSamplesToRead;
	double m_analysisStartTime;
	double m_analysisEndTime;

	long long m_readCycleStartSample;
	long long m_readCycleSamplesToRead;

	unsigned long m_firstMarkedSample;
	unsigned long m_lastMarkedSample;

	bool m_adjustSamplesToRead = true;;

	WaveletParametersSS* m_waveletParams = new WaveletParametersSS;

	int m_patIdx;
	int m_montageNr;
	bool m_useChannsFileBool;
	bool m_useMontage;

	bool m_compareDetectionsWithVisalMarks;
	bool m_MASS_File;

	long m_minutesAnalyzed;
	bool m_verbose;
	bool m_saveTextFile;

	DateAndTime m_fileStartDateAndTime;
};

