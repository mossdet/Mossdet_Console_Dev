#pragma once

#include <string>
#include "DataTypes.h"
#include "HFO_EventsHandler.h"

typedef struct WaveletParametersHFO {
	long lOverlapSamples = 0;
	unsigned int waveletOscillations = 6;
	double dStartFrequency = 0;
	double dEndFrequency = 0;
	double dDeltaFrequency = 0;
	unsigned nrComponents;
};

typedef struct channelAndResection {
	std::string channelName;
	std::string engelOuctome;
	int soz, electrodeType, artifact, resection;
	double rippleRatePerMinute, frRatePerMinute, hfoRatePerMinute, spikeRatePerMinute;
	double maxRippleRatePerMin, maxFRRatePerMin, maxHFO_RatePerMinute, maxSpikeRatePerMinute;
	double totalRipples, totalFRs, totalHFO, totalSpikes;
	double nrMinutes;

	double spikeOrRippleRatePerMinute, spikeOrFrRatePerMinute, spikeOrHfoRatePerMinute;
	double spikeAndRippleRatePerMinute, spikeAndFrRatePerMinute, spikeAndHfoRatePerMinute;
	double spikeAndRipplePlusFRs_RatePerMinute;
	double spindleRippleRatePerMinute, spindleFrRatePerMinute, spindleHfoRatePerMinute;

	double totalSpikeOrRipple, totalSpikeOrFR, totalSpikeOrHFO;
	double totalSpikeAndRipple, totalSpikeAndFR, totalSpikeAndHFO;
	double totalSpikeAndRipplePlusFRs;
	double totalSpindleRipple, totalSpindleFR, totalSpindleHFO;

	/*double auc_rippleRatePerMinute, auc_frRatePerMinute, auc_hfoRatePerMinute;
	double auc_maxRippleRatePerMin, auc_maxFRRatePerMin, auc_maxHFO_RatePerMinute;
	double auc_totalRipples, auc_totalFRs, auc_totalHFO;*/


	void clear() {
		channelName.clear();
		soz = 0;  electrodeType = 0;  artifact = 0;  resection = 0;
		rippleRatePerMinute = 0; frRatePerMinute = 0; hfoRatePerMinute = 0, spikeRatePerMinute = 0;
		maxRippleRatePerMin = 0; maxFRRatePerMin = 0; maxHFO_RatePerMinute = 0, maxSpikeRatePerMinute = 0;
		totalRipples = 0; totalFRs = 0; totalHFO = 0, totalSpikes = 0;
		nrMinutes = 0;

		spikeOrRippleRatePerMinute = 0, spikeOrFrRatePerMinute = 0, spikeOrHfoRatePerMinute = 0;
		spikeAndRippleRatePerMinute = 0, spikeAndFrRatePerMinute = 0, spikeAndHfoRatePerMinute = 0; spikeAndRipplePlusFRs_RatePerMinute = 0;
		spindleRippleRatePerMinute = 0, spindleFrRatePerMinute = 0, spindleHfoRatePerMinute = 0;

		totalSpikeOrRipple = 0, totalSpikeOrFR = 0, totalSpikeOrHFO = 0;
		totalSpikeAndRipple = 0, totalSpikeAndFR = 0, totalSpikeAndHFO = 0;
		totalSpikeAndRipplePlusFRs = 0;
		totalSpindleRipple = 0, totalSpindleFR = 0, totalSpindleHFO = 0;

		/*auc_rippleRatePerMinute = 0; auc_frRatePerMinute = 0; auc_hfoRatePerMinute = 0;
		auc_maxRippleRatePerMin = 0; auc_maxFRRatePerMin = 0; auc_maxHFO_RatePerMinute = 0;
		auc_totalRipples = 0; auc_totalFRs = 0; auc_totalHFO = 0;*/
	}
};

class HFO_Detector
{
public:

	//HFO_Detector(System::String^ strInputFileName, System::String^ strDirName, int patIdx, int montage, std::vector<std::string> *selectedUnipolarChanels, std::vector<std::string> *selectedBipolarChanels, std::string &eoiName, bool useChannsFileBool, System::Double startTime, System::Double endTime, bool compareDetectionsWithVisalMarks);
	HFO_Detector(std::string strInputFileName, std::string detectorPath, std::string strDirName, int patIdx, int montage, std::vector<std::string> *selectedUnipolarChanels,
		std::vector<std::string> *selectedBipolarChanels, std::string &eoiName, bool useChannsFileBool, double startTime, double endTime, bool compareDetectionsWithVisalMarks, double samplingRate, bool verbose, bool saveTextFile);

	long readEEG();
	int characterizeEEG();
	bool detect(EOI_Features &signalFeatures, long firstLocalSampleToRead);

	int getWaveletCausedDelay(void);
	int getFilteringCausedDelay(void);
	bool getFeatureAnalysisAllSignals(matrixStd &rawSignal, EOI_Features &epochFeatures, long rewind, long overlap, long firstSampleToRead, long &samplesToRead, long epochLength, int signalType);
	void getSignalParams(double &samplingRate, double &analysisDuration, int &nrChannels);

	bool getWaveletFrequencyeAnalysisAllSignal(matrixStd &rawSignal, EOI_Features &epochFeatures, long rewind, long overlap, long firstSampleToRead, long &samplesToRead,
		int firstChannel, int lastChannel, long epochLength);

	std::string getPatientNameAndPath();
	double executionTimeLog(int loopNr, long firstSampleToRead, long samplesToRead, long unreadSamples, long totalFileSamples, double durationMs);
	bool saveSelectedEEG_DataToFile(long firstSampleToRead, unsigned rewind, unsigned overlap, matrixStd& signal, matrixStd& notchedSignal, matrixStd& bpSignal, unsigned channelToSave);
	bool plotSelectedEE_Channel(long firstSampleToRead, unsigned rewind, unsigned overlap, matrixStd& signal, matrixStd& notchedSignal, matrixStd& rippleSignal, matrixStd& fastRippleSignal, unsigned channelToSave);



	void readChannelsFromFile(void);
	bool generateEpochCharacterizationRipples(EOI_Features &epochFeatures);
	bool generateEpochCharacterizationFastRipples(EOI_Features &epochFeatures);
	void readAndWriteEOI_VisualAnnotations();
	void getSpikesAndSpikeHFOs(void);
	bool getChannelNumber(std::string markChannelName, unsigned int &channNr);
	bool getNeighboringNumber(long annotChannNr, long &lowerNeighborChNr, long &upperNeighborChNr);
	void generate_EOI_ChannelsFile(void);
	void generate_SpikeChannelsFile(void);
	bool getEpochMarkCoOccurrence(EOI_Features &epochFeatures, int eventType);
	bool areEventsOverlapping(double startTimeA, double endTimeA, double startTimeB, double endTimeB, double minSharedPercent);
	int getHFO_Code(std::string mark);
	void getSpecificStartAndEndTimes();

	bool generateEpochCharacterizationSpikes(EOI_Features &epochFeatures);
	bool generateEpochCharacterizationSpikesHFO(EOI_Features &epochFeatures);

	void writeVectorToFile(int channNr, double time, std::vector<double> frequency, std::vector<double> vectorB);

	void readChannelsAndResection();
	void getChannelsHFO_Statistics(std::string channelName, std::vector<EOI_Event> &readEvents);
	void getROCs(bool localizeResectionOrSOZ);
	void getPercentageROCs(bool localizeResectionOrSOZ);

	void setBiomarkerOccChannels(void);
	void getBiomarkerOcurrenceRates(int channIdx, std::vector<EOI_Event> &readEvents);
	bool generateFileBiomarkersStatistics();	

	void writeInterPatNormParamVals_File(bool localizeResectionOrSOZ, int minute, std::string paramName, std::vector<bool> resctionPerChannel, std::vector<double> paraValsPerChann);
	void writeROC_File(bool localizeResectionOrSOZ, std::string parameterName, int minute, std::vector<double> thVec, std::vector<double> fprVec, std::vector<double> sensVec, std::vector<double> precVec, std::vector<double> mccVector, std::vector<double> kappaVector, double areUnderCurve);
	void writeAUC_File(bool localizeResectionOrSOZ, int minute, std::string parameterName, double areUnderCurve, double bestMCC, double maxMCC_Sens, double maxMCC_Spec, double maxMCC_Prec, double maxMCC_Th, double maxKappa);
	void writeCommonAUC_File(bool localizeResectionOrSOZ, int minute, std::string parameterName, double areUnderCurve, double bestMCC, double maxMCC_Sens, double maxMCC_Spec, double maxMCC_Prec, double maxMCC_Th, double maxKappa);
	void getMaxThAndStepSize_ROC(int paramIdx, double &maxTh, double &stepSize);
	void getChannelParameters();


	void writeDetectionsForSections();
	void readDetectionsForSections();
	void generatePerChannelAnnotations();
	bool writePerChannelDetectionsFile(std::string directory, std::string eventType, std::string channName, std::vector<EOI_Event> &detections);

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

	bool fileExists(const char *fileName) {
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

	std::vector<std::vector<EOI_Event>> m_allMarks;
	std::vector<std::vector<EOI_Event>> m_markedEOI;
	std::vector<std::vector<EOI_Event>> m_markedSpikesHFO;
	std::vector<std::vector<EOI_Event>> m_markedSpikesAlone;
	std::vector<std::vector<EOI_Event>> m_markedAllSpikes;

	double m_samplingRate;
	int m_nrChannels;
	long long m_totalSampleCount;
	long long m_fileAnalysisStartSample;
	long long m_fileAnalysisSamplesToRead;
	double m_analysisStartTime;
	double m_analysisEndTime;

	long long m_readCycleStartSample;
	long long m_readCycleSamplesToRead;

	unsigned long m_firstMarkedSample;
	unsigned long m_lastMarkedSample;

	bool m_adjustSamplesToRead = true;

	WaveletParametersHFO *m_waveletParams = new WaveletParametersHFO;

	int m_patIdx;
	int m_montageNr;
	bool m_useChannsFileBool;
	bool m_useMontage;

	bool m_compareDetectionsWithVisalMarks;

	std::vector<channelAndResection> m_channelsHFO_Statistics;
	std::vector<channelAndResection> m_biomarkerOccPerChann;

	DateAndTime m_fileStartDateAndTime;

	matrixStd ripplesPerChPerMin;
	matrixStd frOccPerMin;
	matrixStd hfoOccPerMin;
	matrixStd iesOccPerMin;
	matrixStd iesRippOccPerMin;
	matrixStd iesFR_OccPerMin;
	matrixStd iesHFO_OccPerMin;

	long m_minutesAnalyzed;

	std::vector<EOI_Event> m_allDetections;

	bool m_verbose;
	bool m_saveTextFile;

};

