#include <chrono>
#include <random>
#include <fstream>
#include <direct.h>
#include <windows.h>
#include <math.h>
#include <map>;
#include <unordered_map>
#include <ctime>
#include <time.h>
#include <iostream>

#include <omp.h>

#include "HFO_Detector.h"
#include "IO_EDFPlus_SignalImporter_c.h"
#include "HFO_EventsHandler.h"
#include "SignalProcessing.h"
#include "CorrelationFeatureSelection.h"
#include "SVM_Detector.h"
#include "Plotter.h"
#include "MatlabFilesRW.h"

#include "Tests.h"

#define PI 3.14159265

#define RPPL_LOW_FREQ			80
#define RPPL_HIGH_FREQ			250
#define FAST_RPPL_LOW_FREQ		250
#define FAST_RPPL_HIGH_FREQ		500

//Epoch labels
#define RAW						1
#define NOTCHED					2
#define RIPPLE					3
#define FAST_RIPPLE				4

#define EPOCH				0.025//s

#define MEASURE_BIOMARKER_VALUE false
#define CALGARY	false

//#define DETECT_ONLY_HFO
#define OCURRENCE_RATE


HFO_Detector::HFO_Detector(std::string strInputFileName, std::string detectorPath, std::string strDirName, int patIdx, int montage, std::vector<std::string> *selectedUnipolarChanels, std::vector<std::string> *selectedBipolarChanels,
	std::string &eoiName, bool useChannsFileBool, double startTime, double endTime, bool compareDetectionsWithVisalMarks, double samplingRate, bool verbose, bool saveTextFile) {

	m_patIdx = patIdx;

	m_strFileName = strInputFileName;
	m_exePath = detectorPath;
	m_outputDirectory = strDirName;

	getPatientNameAndPath();

	m_montageNr = montage;
	m_selectedUnipolarChanels = *selectedUnipolarChanels;
	m_selectedBipolarChanels = *selectedBipolarChanels;
	m_useMontage = montage > 0;

	m_eoiName = eoiName;
	m_useChannsFileBool = useChannsFileBool;
	m_analysisStartTime = startTime;
	m_analysisEndTime = endTime;

	m_waveletParams->dStartFrequency = 20;
	m_waveletParams->dEndFrequency = 520;
	m_waveletParams->dDeltaFrequency = 20;
	m_waveletParams->waveletOscillations = 6;
	//m_exePath = ExePath();

	m_firstMarkedSample = -1;
	m_lastMarkedSample = 0;

	m_compareDetectionsWithVisalMarks = compareDetectionsWithVisalMarks;

	m_minutesAnalyzed = 0;

	m_verbose = verbose;
	m_saveTextFile = saveTextFile;

	double maxNumThreads = omp_get_max_threads();
	maxNumThreads *= 0.75;
	omp_set_num_threads((int)maxNumThreads);
	m_samplingRate = samplingRate;
}


long HFO_Detector::readEEG() {
	
	return m_totalSampleCount;
}

int HFO_Detector::characterizeEEG() {
	if (m_verbose)
		std::cout << "characterizeEEG" << std::endl;
	
	MatlabFilesRW matlabIO_Engine(m_strFileName, m_samplingRate, m_outputDirectory, m_verbose);
	m_nrChannels = 1;
	ContactNames channName;
	channName.contactName = m_patientName;
	channName.contactGlobalIdx = 1;
	m_unipolarLabels.push_back(channName);

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	_mkdir(m_outputDirectory.c_str());
	_mkdir(outputPath.c_str());

	if(m_saveTextFile)
		generate_EOI_ChannelsFile();

	/*******************************************************************************************************************************************************************************************************************/

	long totalFileSamples = matlabIO_Engine.getSampleCount();
	if ((m_analysisStartTime*m_samplingRate) > totalFileSamples)
		m_analysisStartTime = (double)totalFileSamples / m_samplingRate;
	if ((m_analysisEndTime*m_samplingRate) > totalFileSamples)
		m_analysisEndTime = (double)totalFileSamples / m_samplingRate;

	m_fileAnalysisStartSample = (m_analysisStartTime * m_samplingRate) > 0 ? (m_analysisStartTime * m_samplingRate) : 0;
	m_fileAnalysisSamplesToRead = (m_analysisEndTime - m_analysisStartTime) * m_samplingRate;// m_samplingRate.getSampleCount();
	int matlabSC = totalFileSamples;
	if(m_fileAnalysisSamplesToRead > totalFileSamples)
		m_fileAnalysisSamplesToRead = totalFileSamples;
	
	m_totalSampleCount = m_fileAnalysisSamplesToRead;
	//m_fileStartDateAndTime = NULL; m_samplingRate.getDateAndTime();

	long long bigUnreadSamples = m_fileAnalysisSamplesToRead;
	m_readCycleSamplesToRead = m_samplingRate * 60 * 60;		//Detection and data read performed every 10 min , m_fileAnalysisSamplesToRead
	if (m_readCycleSamplesToRead > m_fileAnalysisSamplesToRead)
		m_readCycleSamplesToRead = m_fileAnalysisSamplesToRead;

	m_readCycleStartSample = m_fileAnalysisStartSample;
	long readCycleReadSamples = 0;

	while (bigUnreadSamples > 0) {
		m_reformattedSignalSamples.clear();
		if (m_useMontage) {
			//bool ok = m_samplingRate.ReadReformattedSamples(m_readCycleStartSample, m_readCycleSamplesToRead, m_reformattedSignalSamples);
			//while (!ok);
			;
		}
		else {
			bool ok = matlabIO_Engine.readMatlabSamples(m_readCycleStartSample, m_readCycleSamplesToRead, m_reformattedSignalSamples);
			while (!ok);
		}

		readCycleReadSamples = m_reformattedSignalSamples.size() > 0 ? m_reformattedSignalSamples.back().size() : 0;
		double epochLength = ceil(EPOCH * m_samplingRate);		//20ms Epoch Length
		//std::cout << "readCycleReadSamples" << readCycleReadSamples << std::endl;

		double normLength = 60;	//Normalization every 10 s, 1*60.0
#ifdef OCURRENCE_RATE
		normLength = 60;
#endif
		long processCycleSampleCount = readCycleReadSamples-1;
		long samplesToProcess = normLength * m_samplingRate <= processCycleSampleCount ? normLength * m_samplingRate : processCycleSampleCount;
		long firstSampleToProcess = 0;
		long unprocessedSamples = processCycleSampleCount;
		int loopNr = 0;

		while (unprocessedSamples > 0) {

			m_minutesAnalyzed++;

			//auto start_time = std::chrono::high_resolution_clock::now();

			EOI_Features signalFeatures;
			signalFeatures.setMatrices(m_nrChannels);

			int waveletDelay = getWaveletCausedDelay();// +2 * epochLength;
			int filterDelay = 2 * getFilteringCausedDelay();// +2 * epochLength;
			int delay = waveletDelay >= filterDelay ? 1 * waveletDelay : 1 * filterDelay;

			unsigned rewind = 0;
			if (firstSampleToProcess > delay) {
				rewind = delay;
			}
			unsigned overlap = 0;
			if (firstSampleToProcess + samplesToProcess + delay < processCycleSampleCount) {
				overlap = delay;
			}

			long cycleReadStart = firstSampleToProcess - rewind;
			long cycleReadEnd = firstSampleToProcess + samplesToProcess + overlap;
			long totalReadSymplesCycle = m_reformattedSignalSamples.size();
			//Read cyccle samples
			matrixStd cycleReformattedSamples(m_reformattedSignalSamples.size());
			
#pragma omp parallel for
			for (int ch = 0; ch < m_reformattedSignalSamples.size(); ++ch) {
				for (long s = cycleReadStart; s < cycleReadEnd; ++s) {
					double val = CALGARY ? m_reformattedSignalSamples[ch][s] * -1 : m_reformattedSignalSamples[ch][s];
					cycleReformattedSamples[ch].push_back(val);
					//cycleReformattedSamples[ch].push_back(0);
				}
				//cycleReformattedSamples[ch][(cycleReadEnd - cycleReadStart)/2] = 1;
			}

			matrixStd rawSignal;
			//Filter selected channels
			matrixStd antiOffset_plus_notch, rippleSelctdChs, fastRippleSelectedChs;

			matrixStd matA;
			SignalProcessing::applyHighPassFilterToSignal(cycleReformattedSamples, antiOffset_plus_notch, 1, m_samplingRate);										// Filter offset and DC components
																																									//SignalProcessing::applyBandStopFilterToSignal(matA, antiOffset_plus_notch, 49, 51, m_samplingRate); matA.clear();							//50

																																									// Filter in the frequency band of EOI
			SignalProcessing::rippleBP(cycleReformattedSamples, rippleSelctdChs, RPPL_LOW_FREQ, RPPL_HIGH_FREQ, m_samplingRate);
			SignalProcessing::fastRippleBP(cycleReformattedSamples, fastRippleSelectedChs, FAST_RPPL_LOW_FREQ, FAST_RPPL_HIGH_FREQ, m_samplingRate);

			//saveSelectedEEG_DataToFile(firstSampleToProcess, rewind, overlap, cycleReformattedSamples, antiOffset_plus_notch, rippleSelctdChs, 0);
			//saveSelectedEEG_DataToFile(firstSampleToProcess, rewind, overlap, cycleReformattedSamples, antiOffset_plus_notch, fastRippleSelectedChs, 0);
			//for (int ch = 0; ch < m_reformattedSignalSamples.size(); ++ch) {
			//saveSelectedEEG_DataToFile(firstSampleToProcess, rewind, overlap, cycleReformattedSamples, antiOffset_plus_notch, rippleSelctdChs, ch);
			//plotSelectedEE_Channel(firstSampleToProcess, rewind, overlap, cycleReformattedSamples, antiOffset_plus_notch, rippleSelctdChs, fastRippleSelectedChs, ch);
			//}

			matrixStd waveletPowOutput;
			int firstChann = 0;
			int lastChann = cycleReformattedSamples.size() - 1;
			bool waveletError = getWaveletFrequencyeAnalysisAllSignal(cycleReformattedSamples, signalFeatures, rewind, overlap, firstSampleToProcess, samplesToProcess, firstChann, lastChann, epochLength);

#pragma omp parallel sections
			{
#pragma omp section
				{
					getFeatureAnalysisAllSignals(cycleReformattedSamples, signalFeatures, rewind, overlap, firstSampleToProcess, samplesToProcess, epochLength, RAW);
				}

#pragma omp section
				{
					getFeatureAnalysisAllSignals(antiOffset_plus_notch, signalFeatures, rewind, overlap, firstSampleToProcess, samplesToProcess, epochLength, NOTCHED);
				}
#pragma omp section
				{
					getFeatureAnalysisAllSignals(rippleSelctdChs, signalFeatures, rewind, overlap, firstSampleToProcess, samplesToProcess, epochLength, RIPPLE);
				}
#pragma omp section
				{
					getFeatureAnalysisAllSignals(fastRippleSelectedChs, signalFeatures, rewind, overlap, firstSampleToProcess, samplesToProcess, epochLength, FAST_RIPPLE);
				}
			}

			if (signalFeatures.m_time.empty())
				return -1; 

			bool detectionOK = detect(signalFeatures, firstSampleToProcess);
			if (!detectionOK)
				return -2;

			unprocessedSamples -= samplesToProcess;
			firstSampleToProcess += samplesToProcess;
		}


		bigUnreadSamples -= readCycleReadSamples;
		m_readCycleStartSample += readCycleReadSamples;
	}
	
	/*if (m_verbose)
		std::cout << "writeMatlabSamples" << std::endl;*/
	
	//matlabIO_Engine.writeMatlabSamples(m_allDetections);
	//generatePerChannelAnnotations();

	return 1;
}

bool HFO_Detector::detect(EOI_Features &signalFeatures, long firstLocalSampleToRead) {

	int nrChanns = signalFeatures.features[0]->size();
	int epochSkip = firstLocalSampleToRead == 0 ? 100 : 0;		//skip 250ms, avoid strange feature values stemming from big dc shifts foudn at teh beginning of some signals
	double maxOnsetDiff = 1000;
	double minSharedPercent = 0.75; 0.1; 75;
	bool joinVisualMarks = true;
	bool useWindowConfInterv = false;// (m_patientName.find("SNR0dB") != std::string::npos) || (m_patientName.find("SNR5dB") != std::string::npos) || (m_patientName.find("SNR10dB") != std::string::npos) || (m_patientName.find("SNR15dB") != std::string::npos);

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = m_outputDirectory + "/MOSSDET_Output/" + m_patientName;
	
	_mkdir(patPath.c_str());
	
	for (int channIdx = 0; channIdx < nrChanns; ++channIdx) {
		std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
		int idNumber = m_useMontage ? m_montageLabels[channIdx].montageMOSSDET_Nr : m_unipolarLabels[channIdx].contactGlobalIdx;

		//Fill training epochs
		long nrEpochs = signalFeatures.features[0]->at(channIdx).size();
		std::vector<EOI_Epoch>  inputEpochs(nrEpochs - epochSkip);
		for (int epIdx = epochSkip; epIdx < nrEpochs; ++epIdx) {
			inputEpochs[epIdx - epochSkip].channel = channIdx;
			inputEpochs[epIdx - epochSkip].time = signalFeatures.m_time[epIdx];

			int numFeats = signalFeatures.features.size();
			inputEpochs[epIdx - epochSkip].inputVals.resize(numFeats);

			for (int featIdx = 0; featIdx < numFeats; ++featIdx)
				inputEpochs[epIdx - epochSkip].inputVals[featIdx] = signalFeatures.features[featIdx]->at(channIdx).at(epIdx);
		}

		///AR
		std::vector<unsigned> rippleSelectedFeatures{ 4, 5, 14, 22, 28, 36, 40, 41, 42, 43, 45, 46, 48, 49, 50 };
		std::string rippleFunctionName = m_exePath + "/svm_ripple.dat";
		double rippleTh = 4.6;

		std::vector<unsigned> fastRippleSelectedFeatures{ 8, 9, 14, 22, 28, 36, 54, 55, 56, 57, 64 };
		std::string fastRippleFunctionName = m_exePath + "/svm_fr.dat";
		double fastRippleTh = 1.6;

		std::vector<unsigned> spikeSelectedFeatures{ 12, 14, 15, 23, 24, 26, 28, 29, 31, 37 };
		std::string spikeFunctionName = m_exePath + "/svm_ies.dat";
		double spikeTh = 0.0;

		SVM_Detector svmHFO_Detector("AllHFO", m_patientName, rippleFunctionName, "", rippleSelectedFeatures, m_outputDirectory, m_saveTextFile);
		svmHFO_Detector.setInputEvents(inputEpochs);
		svmHFO_Detector.setAbsoluteSignalLength(m_totalSampleCount / m_samplingRate);
		svmHFO_Detector.setPatientPath(m_patientPath);
		svmHFO_Detector.setSamplingRate(m_samplingRate);
		svmHFO_Detector.setCI(useWindowConfInterv);
		svmHFO_Detector.setFileStartDateAndTime(m_fileStartDateAndTime);

		// Ripple Detection
		std::vector<EOI_Event> readRippleEvents;
		svmHFO_Detector.setEOI_Name("Ripple");
		svmHFO_Detector.setSelectedFeatures(rippleSelectedFeatures);
		svmHFO_Detector.setFunctionName(rippleFunctionName);
		svmHFO_Detector.setProbThreshold(0, rippleTh);
		svmHFO_Detector.clearDetectedEpochs();
		svmHFO_Detector.initEpochVectors(1);
		svmHFO_Detector.detectEpochs(0);
		svmHFO_Detector.getEEG_EventsFromDetectedEpochs(1, channelName);
		svmHFO_Detector.getEvents(readRippleEvents);
		svmHFO_Detector.generateDetectionsAndFeaturesFile(maxOnsetDiff, minSharedPercent);
		svmHFO_Detector.clearEpochVectors();

		// Fast-Ripple Detection
		std::vector<EOI_Event> readFastRippleEvents;
		svmHFO_Detector.setEOI_Name("FastRipple");
		svmHFO_Detector.setSelectedFeatures(fastRippleSelectedFeatures);
		svmHFO_Detector.setFunctionName(fastRippleFunctionName);
		svmHFO_Detector.setProbThreshold(0, fastRippleTh);
		svmHFO_Detector.initEpochVectors(1);
		svmHFO_Detector.detectEpochs(0);
		svmHFO_Detector.getEEG_EventsFromDetectedEpochs(1, channelName);
		svmHFO_Detector.getEvents(readFastRippleEvents);
		svmHFO_Detector.generateDetectionsAndFeaturesFile(maxOnsetDiff, minSharedPercent);
		svmHFO_Detector.clearEpochVectors();

		// Spike Detection
		std::vector<EOI_Event> readSpikeEvents;
		svmHFO_Detector.setEOI_Name("Spike");
		svmHFO_Detector.setSelectedFeatures(spikeSelectedFeatures);
		svmHFO_Detector.setFunctionName(spikeFunctionName);
		svmHFO_Detector.setProbThreshold(0, spikeTh);
		svmHFO_Detector.clearDetectedEpochs();
		svmHFO_Detector.initEpochVectors(1);
		svmHFO_Detector.detectEpochs(0);
		svmHFO_Detector.getEEG_EventsFromDetectedEpochs(1, channelName);
		svmHFO_Detector.getEvents(readSpikeEvents);
		svmHFO_Detector.generateDetectionsAndFeaturesFile(maxOnsetDiff, minSharedPercent);
		svmHFO_Detector.clearEpochVectors();
	}

	return true;
}

bool HFO_Detector::getWaveletFrequencyeAnalysisAllSignal(matrixStd &rawSignal, EOI_Features &epochFeatures, long rewind, long overlap, long firstSampleToProcess, long &samplesToProcess,
	int firstChannel, int lastChannel, long epochLength) {

	double startFrequency = m_waveletParams->dStartFrequency;
	double endFrequency = m_waveletParams->dEndFrequency;
	double deltaFrequency = m_waveletParams->dDeltaFrequency;
	int nrOscillations = m_waveletParams->waveletOscillations;

	long signalLength = rawSignal.back().size();
	long analysisStart = rewind;
	long analysisEnd = signalLength - overlap;

	//get data only from the channels to analyze
	int nrChannsToAnalyze = lastChannel - firstChannel + 1;
	
	matrixStd waveletPowOutput;
	unsigned nrComponents = SignalProcessing::morletWaveletTransform(rawSignal, m_samplingRate, startFrequency, endFrequency, deltaFrequency, nrOscillations, waveletPowOutput);
	SignalProcessing::getMatrixPower(waveletPowOutput);

	unsigned expectedComponents = (endFrequency - startFrequency) / deltaFrequency;
	if (nrComponents != expectedComponents)
		return false;

#pragma omp parallel for
	for (int channIdx = 0; channIdx < nrChannsToAnalyze; ++channIdx) {
		int epochAvgSampleIdx = 0;

		double sumEpochSubHFOPower = 0, sumEpochRpplPower = 0, sumEpochFastRpplPower = 0;
		double avgEpochMaxSubHFOPower = 0, avgEpochMaxRpplPower = 0, avgEpochMaxFastRpplPower = 0;
		double avgEpochMinSubHFOPower = 0, avgEpochMinRpplPower = 0, avgEpochMinFastRpplPower = 0;
		double avgEpochSptrlPeakSubHFO = 0, avgEpochSpctrlPeakRipple = 0, avgEpochSpectralPeakFR = 0;
		double avgEpochMeanMaxRatioSubHFO = 0, avgEpochMeanMaxRatioRipple = 0, avgEpochMeanMaxRatioFR = 0;

		for (long sampIdx = analysisStart; sampIdx < analysisEnd; ++sampIdx) {
			double sumSmplSubHFOPower = 0, sumSmplRpplPower = 0, sumSmplFastRpplPower = 0;
			double maxSmplSubHFOPower = 0, maxSmplRpplPower = 0, maxSmplFastRpplPower = 0;
			double minSmplSubHFOPower = 1e5, minSmplRpplPower = 1e5, minSmplFastRpplPower = 1e5;
			double smplSpectralPeakSubHFO = 0, smplSpectralPeakRipple = 0, smplSpectralPeakFR = 0;
			double smplMeanMaxRatioSubHFO = 0, smplMeanMaxRatioRipple = 0, smplMeanMaxRatioFR = 0;
			std::vector<double> vectorA, vectorB;

			for (unsigned int nc = 0; nc < nrComponents; nc++) {
				int freqBand = startFrequency + deltaFrequency * nc;

				double centerFrequency = startFrequency + nc * deltaFrequency;
				long kernelLenght = centerFrequency > 0 ? nrOscillations / centerFrequency * m_samplingRate : 0;
				double delay = kernelLenght / 2;

				int delayIdx = sampIdx + delay < signalLength ? sampIdx + delay : signalLength - 1;
				double freqPower = waveletPowOutput[channIdx * nrComponents + nc][delayIdx];

				if (!std::isfinite(freqPower) || std::isnan(freqPower) != 0) {
					freqPower = 0;
				}

				if (freqBand >= 0 && freqBand < 79.9) {
					sumSmplSubHFOPower += freqPower;
					if (freqPower > maxSmplSubHFOPower) {
						maxSmplSubHFOPower = freqPower;
						smplSpectralPeakSubHFO = freqBand;
					}
					if (freqPower < minSmplSubHFOPower)
						minSmplSubHFOPower = freqPower;
				}
				else if (freqBand >= 79.9 && freqBand <= 251.1) {
					sumSmplRpplPower += freqPower;
					if (freqPower > maxSmplRpplPower) {
						maxSmplRpplPower = freqPower;
						smplSpectralPeakRipple = freqBand;
					}
					if (freqPower < minSmplRpplPower)
						minSmplRpplPower = freqPower;
				}
				else if (freqBand >= 251.1) {
					sumSmplFastRpplPower += freqPower;
					if (freqPower > maxSmplFastRpplPower) {
						maxSmplFastRpplPower = freqPower;
						smplSpectralPeakFR = freqBand;
					}
					if (freqPower < minSmplFastRpplPower)
						minSmplFastRpplPower = freqPower;
				}

				vectorA.push_back(freqBand);
				vectorB.push_back(freqPower);

			}

			//writeVectorToFile(channIdx, ((double)(m_fileAnalysisStartSample + firstSampleToProcess + sampIdx - (rewind)+1) / m_samplingRate), vectorA, vectorB);
			
			smplMeanMaxRatioSubHFO = (sumSmplSubHFOPower / (double)nrComponents) / maxSmplSubHFOPower;	//max - diff
			smplMeanMaxRatioRipple = (sumSmplRpplPower / (double)nrComponents) / maxSmplRpplPower;		//max - diff
			smplMeanMaxRatioFR = (sumSmplFastRpplPower / (double)nrComponents) / maxSmplFastRpplPower;	//max - diff

			sumEpochSubHFOPower += sumSmplSubHFOPower; sumEpochRpplPower += sumSmplRpplPower; sumEpochFastRpplPower += sumSmplFastRpplPower;
			avgEpochMaxSubHFOPower += maxSmplSubHFOPower; avgEpochMaxRpplPower += maxSmplRpplPower; avgEpochMaxFastRpplPower += maxSmplFastRpplPower;
			avgEpochMinSubHFOPower += minSmplSubHFOPower; avgEpochMinRpplPower += minSmplRpplPower; avgEpochMinFastRpplPower += minSmplFastRpplPower;
			avgEpochSptrlPeakSubHFO += smplSpectralPeakSubHFO; avgEpochSpctrlPeakRipple += smplSpectralPeakRipple; avgEpochSpectralPeakFR += smplSpectralPeakFR;
			avgEpochMeanMaxRatioSubHFO += smplMeanMaxRatioSubHFO; avgEpochMeanMaxRatioRipple += smplMeanMaxRatioRipple; avgEpochMeanMaxRatioFR += smplMeanMaxRatioFR;

			epochAvgSampleIdx++;
			if (epochAvgSampleIdx == epochLength) {
				sumEpochSubHFOPower /= epochLength; sumEpochRpplPower /= epochLength; sumEpochFastRpplPower /= epochLength;
				avgEpochMaxSubHFOPower /= epochLength; avgEpochMaxRpplPower /= epochLength; avgEpochMaxFastRpplPower /= epochLength;
				avgEpochMinSubHFOPower /= epochLength; avgEpochMinRpplPower /= epochLength; avgEpochMinFastRpplPower /= epochLength;
				avgEpochSptrlPeakSubHFO /= epochLength; avgEpochSpctrlPeakRipple /= epochLength; avgEpochSpectralPeakFR /= epochLength;
				avgEpochMeanMaxRatioSubHFO /= epochLength; avgEpochMeanMaxRatioRipple /= epochLength; avgEpochMeanMaxRatioFR /= epochLength;

				epochFeatures.m_sumEpochSubHFOPower[firstChannel + channIdx].push_back(sumEpochSubHFOPower);
				epochFeatures.m_avgEpochMaxSubHFOPower[firstChannel + channIdx].push_back(avgEpochMaxSubHFOPower);
				epochFeatures.m_avgEpochSptrlPeakSubHFO[firstChannel + channIdx].push_back(avgEpochSptrlPeakSubHFO);
				epochFeatures.m_avgEpochMinMaxRatioSubHFO[firstChannel + channIdx].push_back(avgEpochMinSubHFOPower / avgEpochMaxSubHFOPower);

				epochFeatures.m_sumEpochRpplPower[firstChannel + channIdx].push_back(sumEpochRpplPower);
				epochFeatures.m_avgEpochMaxRpplPower[firstChannel + channIdx].push_back(avgEpochMaxRpplPower);
				epochFeatures.m_avgEpochSpctrlPeakRipple[firstChannel + channIdx].push_back(avgEpochSpctrlPeakRipple);
				epochFeatures.m_avgEpochMinMaxRatioRipple[firstChannel + channIdx].push_back(avgEpochMinRpplPower / avgEpochMaxRpplPower);

				epochFeatures.m_sumEpochFastRpplPower[firstChannel + channIdx].push_back(sumEpochFastRpplPower);
				epochFeatures.m_avgEpochMaxFastRpplPower[firstChannel + channIdx].push_back(avgEpochMaxFastRpplPower);
				epochFeatures.m_avgEpochSpectralPeakFR[firstChannel + channIdx].push_back(avgEpochSpectralPeakFR);
				epochFeatures.m_avgEpochMinMaxRatioFR[firstChannel + channIdx].push_back(avgEpochMinFastRpplPower / avgEpochMaxFastRpplPower);
				
				double time = ((double)(m_readCycleStartSample + firstSampleToProcess + sampIdx - (rewind)+1) / m_samplingRate);
				if (firstChannel + channIdx == 0) {
					epochFeatures.m_time.push_back(time);
				}
				epochFeatures.m_EOI_Mark[channIdx].push_back(0);

				epochAvgSampleIdx = 0;
				sumEpochSubHFOPower = sumEpochRpplPower = sumEpochFastRpplPower = 0;
				avgEpochMaxSubHFOPower = avgEpochMaxRpplPower = avgEpochMaxFastRpplPower = 0;
				avgEpochMinSubHFOPower = 0, avgEpochMinRpplPower = 0, avgEpochMinFastRpplPower = 0;
				avgEpochSptrlPeakSubHFO = avgEpochSpctrlPeakRipple = avgEpochSpectralPeakFR = 0;
				avgEpochMeanMaxRatioSubHFO = avgEpochMeanMaxRatioRipple = avgEpochMeanMaxRatioFR = 0;

				if (sampIdx + 1 < analysisEnd) {				//epochs share 50%
					sampIdx = sampIdx - epochLength / 2;
				}
				if (sampIdx + epochLength >= analysisEnd) {		//go back if not enough sampels for one full epoch
					sampIdx = analysisEnd - epochLength;
				}

			}
		}
	}

	return true;
}

bool HFO_Detector::getFeatureAnalysisAllSignals(matrixStd &signal, EOI_Features &epochFeatures, long rewind, long overlap, long firstSampleToProcess, long &samplesToProcess, long epochLength, int signalType) {

	long signalLength = signal[0].size();
	int nrHFO_Channs = m_nrChannels;
	long analysisStart = rewind;
	long analysisEnd = signalLength - overlap;

//#pragma omp parallel for
	for (int channIdx = 0; channIdx < nrHFO_Channs; ++channIdx) {

		int epochAvgSampleIdx = 0;

		//Amplitude Features
		double maxAmpl = 0, amplVariance = 0, meanAmplitude = 0, lineLength = 0, teagerEnergy = 0;
		//Waveform Features
		double zeroCrossingsNr = 0, peaksRateNr = 0;
		// Hjorth Features
		double mobility = 0, complexity = 0;
		// Symmetry
		double symmetry = 0, assymetry = 0;
		double power = 0;

		double binarization = 0;

		//Helpers
		SignalProcessing varianceCalc;
		if (!varianceCalc.setVarianceData(signal[channIdx]))
			bool stop = true;

		double max = -1e9, min = 1e9;

		//Mean amplitude needds to be obtained one epoch in advance in order to calculate the zero crossings
		int meanAmpBfferSize = epochLength;
		for (long i = analysisStart; i < analysisStart + meanAmpBfferSize; ++i) {
			meanAmplitude += signal[channIdx][i];
		}
		meanAmplitude /= meanAmpBfferSize;

		std::vector<double> autoCorrSignal;
		for (long sampIdx = analysisStart; sampIdx < analysisEnd; ++sampIdx) {

			double val = signal[channIdx][sampIdx];

			if (sampIdx >= analysisStart + meanAmpBfferSize && sampIdx + 1 < signalLength) {
				meanAmplitude *= meanAmpBfferSize;
				meanAmplitude = (meanAmplitude + val - signal[channIdx][sampIdx - meanAmpBfferSize])/ meanAmpBfferSize;
			}


			//Max amplitude
			if (val > max)
				max = val;
			if (val < min)
				min = val;

			//Line Length
			if (sampIdx > 0) {
				lineLength += std::abs(signal[channIdx][sampIdx] - signal[channIdx][sampIdx - 1]);
			}

			if (sampIdx > 0 && sampIdx < signalLength - 1) {																			//Teager energy
				teagerEnergy += (val * val) - (signal[channIdx][sampIdx - 1] * signal[channIdx][sampIdx + 1]);
			}

			// Zero crossings
			if (sampIdx > 0) {
				if (signal[channIdx][sampIdx - 1] > meanAmplitude && val <= meanAmplitude) {
					zeroCrossingsNr++;
				}
				else if (signal[channIdx][sampIdx - 1] < meanAmplitude && val >= meanAmplitude) {
					zeroCrossingsNr++;
				}
			}

			/*if (sampIdx > 2 && sampIdx < signalLength - 2) {
				double before = (signal[channIdx][sampIdx - 2] + signal[channIdx][sampIdx - 1] + val) / 3;
				double center = (signal[channIdx][sampIdx - 1] + val + signal[channIdx][sampIdx + 1]) / 3;
				double after = (signal[channIdx][sampIdx + 2] + signal[channIdx][sampIdx + 1] + val) / 3;*/
			double center = 0, before = 0, after = 0;

			/*if (sampIdx >=2 && sampIdx < signalLength - 2) {
				before = (signal[channIdx][sampIdx - 2] + signal[channIdx][sampIdx - 1] + val) / 3;
				center = (signal[channIdx][sampIdx - 1] + val + signal[channIdx][sampIdx + 1]) / 3;
				after = (signal[channIdx][sampIdx + 2] + signal[channIdx][sampIdx + 1] + val) / 3;
			}
			else*/ 
				if (sampIdx >= 1 && sampIdx < signalLength - 1) {
				/*before = (signal[channIdx][sampIdx - 1] + val) / 2;
				center = (signal[channIdx][sampIdx - 1] + val + signal[channIdx][sampIdx + 1]) / 3;
				after = (signal[channIdx][sampIdx + 1] + val) / 2;*/

				before = signal[channIdx][sampIdx - 1];
				center = val;
				after = signal[channIdx][sampIdx + 1];
			}
			if (compareDoublesLarger(center, before) && compareDoublesLarger(center, after) || compareDoublesSmaller(center, before) && compareDoublesSmaller(center, after)) {
				peaksRateNr++;
			}

			power += val * val;

			//Binarization
			if (sampIdx > 0) {
				binarization += (std::abs(signal[channIdx][sampIdx]) >= abs(signal[channIdx][sampIdx - 1]));
			}

			/****************************************************************************************/
			epochAvgSampleIdx++;
			if (epochAvgSampleIdx == epochLength) {

				long epochSampStart = sampIdx - (epochLength - 1);
				long epochSampleEnd = sampIdx;

				maxAmpl = (max - min);
				amplVariance = varianceCalc.runVariance(epochSampStart, epochLength);
				lineLength;
				zeroCrossingsNr;
				teagerEnergy;
				peaksRateNr;
				binarization /= epochLength;
				binarization = 0;

				// Hjorth Features
				std::vector<double> firstDerivative, secondDerivative, thirdDerivative;
				SignalProcessing::firstDerivative(signal[channIdx], firstDerivative, epochSampStart, epochLength, m_samplingRate);
				SignalProcessing::firstDerivative(firstDerivative, secondDerivative, 0, epochLength, m_samplingRate);
				SignalProcessing::firstDerivative(secondDerivative, thirdDerivative, 0, epochLength, m_samplingRate);
				// Mobility
				SignalProcessing varianceFD;
				if (!varianceFD.setVarianceData(firstDerivative))
					bool stop = true;
				double varValFD = varianceFD.runVariance(0, epochLength);
				mobility = sqrt(varValFD / amplVariance);
				// Complexity
				SignalProcessing varianceSD;
				if (!varianceSD.setVarianceData(secondDerivative))
					bool stop = true;
				double varValSD = varianceSD.runVariance(0, epochLength);
				double mobilityFD = varValSD / varValFD;
				complexity = mobilityFD / mobility;

				//Symmetry and Asymmetry
				double center = (epochSampStart + epochLength / 2);
				double numSymm = 0, nummAss = 0, maxSymm = -1 * 1000, maxAsymm = -1 * 1000;
				long i = center, j = center, n = 0;
				while (i < epochSampleEnd && j >= epochSampStart) {
					double symmetric = 0.5 * (signal[channIdx][i] + signal[channIdx][j]);
					double assymetric = 0.5 * (signal[channIdx][i] - signal[channIdx][j]);
					if (abs(symmetric) > maxSymm)
						maxSymm = abs(symmetric);
					if (abs(assymetric) > maxAsymm)
						maxAsymm = abs(assymetric);
					numSymm += symmetric * symmetric;
					nummAss += assymetric * assymetric;
					i++; j--; n++;
				}
				symmetry = (maxSymm*maxSymm * (2 * n + 1)) > 0 ? numSymm / (maxSymm*maxSymm * (2 * n + 1)) : 0;
				assymetry = (maxAsymm*maxAsymm * (2 * n + 1)) > 0 ? nummAss / (maxAsymm*maxAsymm * (2 * n + 1)) : 0;

				double sumFirstDerivative = 0, sumSecondDerivative = 0, sumThirdDerivative = 0;
				int nrDerivSamples = firstDerivative.size();
				for (int derivIdx = 0; derivIdx < nrDerivSamples; derivIdx++) {
					sumFirstDerivative += abs(firstDerivative[derivIdx]);
					sumSecondDerivative += abs(secondDerivative[derivIdx]);
					sumThirdDerivative += abs(thirdDerivative[derivIdx]);
				}
				sumFirstDerivative;
				sumSecondDerivative;
				sumThirdDerivative;

				if (signalType == RAW) {
					epochFeatures.m_maxAmplRaw[channIdx].push_back(maxAmpl);
					epochFeatures.m_amplVarianceRaw[channIdx].push_back(amplVariance);
					epochFeatures.m_lineLengthRaw[channIdx].push_back(lineLength);
					epochFeatures.m_teagerEnergyRaw[channIdx].push_back(teagerEnergy);
					epochFeatures.m_zeroCrossingsNrRaw[channIdx].push_back(zeroCrossingsNr);
					epochFeatures.m_peaksRateNrRaw[channIdx].push_back(peaksRateNr);
					epochFeatures.m_mobilityRaw[channIdx].push_back(mobility);
					epochFeatures.m_complexityRaw[channIdx].push_back(complexity);
					epochFeatures.m_symmetryRaw[channIdx].push_back(symmetry);
					epochFeatures.m_asymmetryRaw[channIdx].push_back(assymetry);
					epochFeatures.m_powerRaw[channIdx].push_back(power);
					epochFeatures.m_sumFirstDerivRaw[channIdx].push_back(sumFirstDerivative);
					epochFeatures.m_sumSecondDerivRaw[channIdx].push_back(sumSecondDerivative);
					epochFeatures.m_sumThirdDerivRaw[channIdx].push_back(sumThirdDerivative);
				}
				else if (signalType == NOTCHED) {
					epochFeatures.m_maxAmplNotchedDC[channIdx].push_back(maxAmpl);
					epochFeatures.m_amplVarianceNotchedDC[channIdx].push_back(amplVariance);
					epochFeatures.m_lineLengthNotchedDC[channIdx].push_back(lineLength);
					epochFeatures.m_teagerEnergyNotchedDC[channIdx].push_back(teagerEnergy);
					epochFeatures.m_zeroCrossingsNrNotchedDC[channIdx].push_back(zeroCrossingsNr);
					epochFeatures.m_peaksRateNrNotchedDC[channIdx].push_back(peaksRateNr);
					epochFeatures.m_mobilityNotchedDC[channIdx].push_back(mobility);
					epochFeatures.m_complexityNotchedDC[channIdx].push_back(complexity);
					epochFeatures.m_symmetryNotchedDC[channIdx].push_back(symmetry);
					epochFeatures.m_asymmetryNotchedDC[channIdx].push_back(assymetry);
					epochFeatures.m_powerNotchedDC[channIdx].push_back(power);
					epochFeatures.m_sumFirstDerivNotchedDC[channIdx].push_back(sumFirstDerivative);
					epochFeatures.m_sumSecondDerivNotchedDC[channIdx].push_back(sumSecondDerivative);
					epochFeatures.m_sumThirdDerivNotchedDC[channIdx].push_back(sumThirdDerivative);

				}
				else if (signalType == RIPPLE) {
					epochFeatures.m_maxAmplRppl[channIdx].push_back(maxAmpl);
					epochFeatures.m_amplVarianceRppl[channIdx].push_back(amplVariance);
					epochFeatures.m_lineLengthRppl[channIdx].push_back(lineLength);
					epochFeatures.m_teagerEnergyRppl[channIdx].push_back(teagerEnergy);
					epochFeatures.m_zeroCrossingsNrRppl[channIdx].push_back(zeroCrossingsNr);
					epochFeatures.m_peaksRateNrRppl[channIdx].push_back(peaksRateNr);
					epochFeatures.m_mobilityRppl[channIdx].push_back(mobility);
					epochFeatures.m_complexityRppl[channIdx].push_back(complexity);
					epochFeatures.m_symmetryRppl[channIdx].push_back(symmetry);
					epochFeatures.m_asymmetryRppl[channIdx].push_back(assymetry);
					epochFeatures.m_powerRppl[channIdx].push_back(power);
					epochFeatures.m_sumFirstDerivRppl[channIdx].push_back(sumFirstDerivative);
					epochFeatures.m_sumSecondDerivRppl[channIdx].push_back(sumSecondDerivative);
					epochFeatures.m_sumThirdDerivRppl[channIdx].push_back(sumThirdDerivative);

				}
				else if (signalType == FAST_RIPPLE) {
					epochFeatures.m_maxAmplFastRppl[channIdx].push_back(maxAmpl);
					epochFeatures.m_amplVarianceFastRppl[channIdx].push_back(amplVariance);
					epochFeatures.m_lineLengthFastRppl[channIdx].push_back(lineLength);
					epochFeatures.m_teagerEnergyFastRppl[channIdx].push_back(teagerEnergy);
					epochFeatures.m_zeroCrossingsNrFastRppl[channIdx].push_back(zeroCrossingsNr);
					epochFeatures.m_peaksRateNrFastRppl[channIdx].push_back(peaksRateNr);
					epochFeatures.m_mobilityFastRppl[channIdx].push_back(mobility);
					epochFeatures.m_complexityFastRppl[channIdx].push_back(complexity);
					epochFeatures.m_symmetryFastRppl[channIdx].push_back(symmetry);
					epochFeatures.m_asymmetryFastRppl[channIdx].push_back(assymetry);
					epochFeatures.m_powerFastRppl[channIdx].push_back(power);
					epochFeatures.m_sumFirstDerivFastRppl[channIdx].push_back(sumFirstDerivative);
					epochFeatures.m_sumSecondDerivFastRppl[channIdx].push_back(sumSecondDerivative);
					epochFeatures.m_sumThirdDerivFastRppl[channIdx].push_back(sumThirdDerivative);

				}

				// clear feature variables
				maxAmpl = 0; amplVariance = 0; lineLength = 0; teagerEnergy = 0;
				zeroCrossingsNr = 0; peaksRateNr = 0;
				mobility = 0; complexity = 0;
				symmetry = 0; assymetry = 0;
				//autoCorrSignal.clear();
				power = 0;
				binarization = 0;

				//clear helpers
				max = -1e9; min = 1e9;

				epochAvgSampleIdx = 0;

				if (sampIdx + 1 < analysisEnd) {				//epochs share 50%
					sampIdx = sampIdx - epochLength / 2;
				}
				if (sampIdx + epochLength >= analysisEnd) {		//go back if not enough sampels for one full epoch
					sampIdx = analysisEnd - epochLength;
				}
			}
		}
	}

	return true;
}

void HFO_Detector::generate_EOI_ChannelsFile(void) {
	std::ofstream eoiChannels;

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName;
	std::string filename = patPath + "/" + m_patientName + "_ChannelsFile.txt";
	//std::string filename = patPath + "/" + m_patientName + "_AllChannels.txt";
	//std::string filename = patPath + "/" + m_patientName + "_EOI_ChannsFile.txt";

	_mkdir(outputPath.c_str());
	_mkdir(patPath.c_str());
	remove(filename.c_str());

	eoiChannels.open(filename);
	eoiChannels << "Channel\t" << "Label\n";

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
		channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
		eoiChannels << channIdx << "\t" << channelName << "\n";
	}
	eoiChannels.close();
}

void HFO_Detector::generate_SpikeChannelsFile(void) {
	m_markedAllSpikes;
	std::ofstream spikeChannels;

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName;
	std::string filename = patPath + "/" + m_patientName + "_SpikeChannels.txt";
	_mkdir(outputPath.c_str());
	_mkdir(patPath.c_str());

	if (!fileExists(filename.c_str())) {
		spikeChannels.open(filename);
		spikeChannels << "Channel\t" << "Label\n";
	}
	else {
		spikeChannels.open(filename, std::ios::app);    // open file for appending
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
		channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());

		if (!m_markedAllSpikes[channIdx].empty()) {
			spikeChannels << channIdx << "\t" << channelName << "\n";
		}


	}

	spikeChannels.close();
}

int HFO_Detector::getWaveletCausedDelay(void) {

	int maxDelay = 0;

	double startFrequency = m_waveletParams->dStartFrequency;
	double endFrequency = m_waveletParams->dEndFrequency;
	double deltaFrequency = m_waveletParams->dDeltaFrequency;
	int nrOscillations = m_waveletParams->waveletOscillations;

	int nrComponents = (endFrequency - startFrequency) / deltaFrequency;
	for (int componentIdx = 0; componentIdx < nrComponents; componentIdx++) {
		double centerFrequency = startFrequency + componentIdx * deltaFrequency;
		long kernelLenght = centerFrequency > 0 ? (nrOscillations / centerFrequency) * m_samplingRate : 0;
		double delay = kernelLenght / 2;
		if (delay > maxDelay)
			maxDelay = delay;
	}

	return maxDelay;
}

int HFO_Detector::getFilteringCausedDelay(void) {

	return SignalProcessing::getOrder(); //we filter two times so teh delay is equal to the filter order
}

std::string HFO_Detector::getPatientNameAndPath() {
	for (int i = m_strFileName.length() - 1; i >= 0; i--) {
		if (m_strFileName[i] == 92) {
			m_patientName = m_strFileName.substr(i + 1, m_strFileName.length() - i);
			m_patientPath = m_strFileName.substr(0, i + 1);
			break;
		}
	}

	m_patientName.pop_back(); m_patientName.pop_back(); m_patientName.pop_back(); m_patientName.pop_back();

	return m_patientName;
}

double HFO_Detector::executionTimeLog(int loopNr, long firstSampleToProcess, long samplesToProcess, long unprocessedSamples, long totalFileSamples, double durationMs) {
	std::ofstream characterizationTimeLog;
	
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName;
	std::string filename = patPath + "/" + m_patientName + "_ExecutionTime.txt";
	_mkdir(outputPath.c_str());
	_mkdir(patPath.c_str());

	if (!fileExists(filename.c_str())) {
		characterizationTimeLog.open(filename);
		characterizationTimeLog << "LoopNr\t" << "firstSampleToProcess\t" << "SamplesToRead\t" << "UnreadSamples\t" << "TotalFileSamples\t" << "Duration (ms) \n";
	}
	else {
		characterizationTimeLog.open(filename, std::ios::app);    // open file for appending
	}

	characterizationTimeLog << loopNr << "\t" << firstSampleToProcess << "\t" << samplesToProcess << "\t" << unprocessedSamples << "\t" << totalFileSamples << "\t" << durationMs << "\n";
	characterizationTimeLog.close();

	return durationMs;
}

bool HFO_Detector::saveSelectedEEG_DataToFile(long firstSampleToProcess, unsigned rewind, unsigned overlap, matrixStd& signal,matrixStd& notchedSignal, matrixStd& bpSignal, unsigned channelToSave) {

	long nrSamples = signal[0].size();
	int channIdx = channelToSave;// m_EOI_Channels->at(channelToSave);
	std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;

	int firOrder = SignalProcessing::getOrder();
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName;
	std::string filename = patPath + "/" + m_patientName + "_" + channelName + "_Fir" + std::to_string(firOrder) + "_EEG_Data.dat";
	_mkdir(outputPath.c_str());
	_mkdir(patPath.c_str());
	//remove(filename.c_str());

	std::ofstream eegDataFile;
	if (!fileExists(filename.c_str())) {
		eegDataFile.open(filename);
		//Header
		eegDataFile << "Time\t" << "Raw\t" << "Notched+High-Pass\t" << "Band-Pass\t" << "BP_TeagerEnergy\t" << "firstDerivative\t" << "secondDerivative\t" << "thirdDerivative\n";
	}
	else {
		eegDataFile.open(filename, std::ios::app);    // open file for appending
	}

	std::vector<double> firstDerivative, secondDerivative, thirdDerivative;
	SignalProcessing::firstDerivative(signal[channIdx], firstDerivative, 0, nrSamples, m_samplingRate);
	SignalProcessing::firstDerivative(firstDerivative, secondDerivative, 0, nrSamples, m_samplingRate);
	SignalProcessing::firstDerivative(secondDerivative, thirdDerivative, 0, nrSamples, m_samplingRate);

	int timeZeroCtr = 0;
	double reverseFactor = 1.0;
	for (long sampIdx = rewind; sampIdx < nrSamples - overlap; sampIdx++) {
		double time = (double)(m_fileAnalysisStartSample + firstSampleToProcess + sampIdx - rewind) / m_samplingRate;
		double rawVal = signal[channelToSave][sampIdx] * reverseFactor;
		double notchedVal = notchedSignal[channelToSave][sampIdx] * reverseFactor;
		double bpValPrev = sampIdx > 0 ? bpSignal[channelToSave][sampIdx - 1] * reverseFactor : bpSignal[channelToSave][sampIdx] * reverseFactor;
		double bpVal = bpSignal[channelToSave][sampIdx] * reverseFactor;
		double bpValPost = sampIdx+1 < nrSamples - overlap ? bpSignal[channelToSave][sampIdx+1] * reverseFactor : bpSignal[channelToSave][sampIdx] * reverseFactor;

		double bpTeagerEnergy = bpVal * bpVal - bpValPrev * bpValPost;
		eegDataFile << time << "\t" << rawVal << "\t" << notchedVal << "\t" << bpVal << "\t" << bpTeagerEnergy << "\t" << firstDerivative[sampIdx] << "\t" << secondDerivative[sampIdx] << "\t" << thirdDerivative[sampIdx] << "\n";
		if (time == 0.0)
			timeZeroCtr++;
	}

	eegDataFile.close();

	return timeZeroCtr == 1;
}

bool HFO_Detector::plotSelectedEE_Channel(long firstSampleToProcess, unsigned rewind, unsigned overlap, matrixStd& signal, matrixStd& notchedSignal, matrixStd& rippleSignal, matrixStd& fastRippleSignal, unsigned channelToSave) {

	long signalLength = signal[0].size();
	long analysisStart = rewind;
	long analysisEnd = signalLength - overlap;
	int channIdx = channelToSave;// m_EOI_Channels->at(channelToSave);
	std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
	std::string patientPath = "FilterCharacterization/" + m_patientName;
	_mkdir("FilterCharacterization/");
	_mkdir(patientPath.c_str());
	
	int timeZeroCtr = 0;
	std::vector<double> time, raw, highpassed, ripple, fastRipple;
	double reverseFactor = 1.0;
	//for (long sampIdx = rewind; sampIdx < nrSamples - overlap; sampIdx++) {
	for (long sampIdx = analysisStart; sampIdx < analysisEnd; ++sampIdx) {
		time.push_back((double)(m_fileAnalysisStartSample + firstSampleToProcess + sampIdx - rewind) / m_samplingRate);
		raw.push_back(signal[channelToSave][sampIdx] * reverseFactor);
		highpassed.push_back(notchedSignal[channelToSave][sampIdx] * reverseFactor);
		ripple.push_back(rippleSignal[channelToSave][sampIdx] * reverseFactor);
		fastRipple.push_back(fastRippleSignal[channelToSave][sampIdx] * reverseFactor);
	}
	
	//std::string wdwTitleBase = patientPath + "/Fir" + std::to_string(SignalProcessing::getOrder());
	std::string wdwTitleBase = patientPath + "/" + channelName;

	std::string wdwTitle = wdwTitleBase + "_Raw_"+ std::to_string(time.back());
	Plotter::plotXY(time, raw, wdwTitle, 0, 0);

	wdwTitle = wdwTitleBase + "_HP_" + std::to_string(time.back());
	Plotter::plotXY(time, highpassed, wdwTitle, 0, 0);

	wdwTitle = wdwTitleBase + "_Ripple_" + std::to_string(time.back());
	Plotter::plotXY(time, ripple, wdwTitle, 0, 0);

	wdwTitle = wdwTitleBase + "_FastRipple_" + std::to_string(time.back());
	Plotter::plotXY(time, fastRipple, wdwTitle, 0, 0);

	return true;
}


bool HFO_Detector::generateEpochCharacterizationRipples(EOI_Features &epochFeatures) {
	
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patOutPath = outputPath + "/" + m_patientName;
	std::string patTrainingFilesPath = outputPath + "/" + m_patientName + "/TrainingFiles";

	_mkdir(outputPath.c_str());
	_mkdir(patOutPath.c_str());
	_mkdir(patTrainingFilesPath.c_str());

	epochFeatures.normalizeFeatures();/// Just for the paper

	std::vector<bool> chnnWithMarks(m_nrChannels);
	bool ripple = true;
	bool fastRipple = false;
	if (ripple) {
		getEpochMarkCoOccurrence(epochFeatures, RIPPLE_LABEL);
		std::string nrPositivesFilename = patOutPath + "/" + m_patientName + "_RipplesPerCh.txt";
		remove(nrPositivesFilename.c_str());
		std::ofstream nrPositivesFile;
		nrPositivesFile.open(nrPositivesFilename);
		nrPositivesFile << "nrMarks" << "\t" << "channelName" << "\n";
		long nrTotalPositives = 0;

	//#pragma omp parallel for reduction(+ : nrTotalPositives)
		for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
			int nrMarks = 0;
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double mark = epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx];
				if (mark > 0) {
					nrMarks++;
				}
			}
			if (nrMarks > 0) {
				nrTotalPositives += nrMarks;
				chnnWithMarks[channIdx] = true;
				std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
				channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
				nrPositivesFile << nrMarks << "\t" << channelName << "\n";
			}

		}
		nrPositivesFile << nrTotalPositives << "\t" << "TotalPositives" << "\n";
		nrPositivesFile.close();
	}

#pragma omp parallel for
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		if (chnnWithMarks[channIdx]) {
			std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
			channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
			int idNumber = m_useMontage ? m_montageLabels[channIdx].montageMOSSDET_Nr : m_unipolarLabels[channIdx].contactGlobalIdx;
			//std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization_Ripple.txt";
			std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization.txt";
			filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());

			std::ofstream trainingFile;
			if (!fileExists(filename.c_str())) {
				trainingFile.open(filename);
				//Header
				trainingFile << "Data \t" << "CH \t" << "Time \t";
				for (unsigned i = 0; i < epochFeatures.featureNames.size(); ++i) {
					trainingFile << epochFeatures.featureNames[i] << "\t";
				}
				trainingFile << "Output \n";
			}
			else {
				trainingFile.open(filename, std::ios::app);    // open file for appending
			}

			trainingFile.precision(32);
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			//Generate Excel readable File
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double sampleTime = epochFeatures.m_time[sampleIdx];
				unsigned currSample = epochFeatures.m_time[sampleIdx] * m_samplingRate;

				trainingFile << "d:\t" << channelName << "\t" << epochFeatures.m_time[sampleIdx] << "\t";

				unsigned nrFeats = epochFeatures.features.size();
				for (unsigned featIdx = 0; featIdx < nrFeats; ++featIdx) {
					unsigned featLength = epochFeatures.features[featIdx]->at(channIdx).size();
					if (featLength > 1 && sampleIdx < featLength) {
						trainingFile << epochFeatures.features[featIdx]->at(channIdx)[sampleIdx] << "\t";
					}
				}

				trainingFile << epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx] << "\n";
			}
			trainingFile.close();
		}
	}

	epochFeatures.normalizeFeatures();
	std::vector<int> featsToExclude;

	// bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes
	CorrelationFeatureSelection featureSelectorRipples(m_nrChannels, patOutPath, m_patientName);
	for (int i = 0; i <= 3; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 7; i <= 11; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 23; i <= 25; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 37; i <= 39; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 51; i <= 67; ++i) {
		featsToExclude.push_back(i);
	}
	featureSelectorRipples.setFeatsToExclude(featsToExclude);
	featureSelectorRipples.getEOI_MarkFeatureCorrelation(true, false, false, false, false, false, epochFeatures);	//Ripples

	return true;
}

bool HFO_Detector::generateEpochCharacterizationFastRipples(EOI_Features &epochFeatures) {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patOutPath = outputPath + "/" + m_patientName;
	std::string patTrainingFilesPath = outputPath + "/" + m_patientName + "/TrainingFiles";

	_mkdir(outputPath.c_str());
	_mkdir(patOutPath.c_str());
	_mkdir(patTrainingFilesPath.c_str());


	std::vector<bool> chnnWithMarks(m_nrChannels);
	bool ripple = true;
	bool fastRipple = false;
	if (ripple) {
		getEpochMarkCoOccurrence(epochFeatures, FAST_RIPPLE_LABEL);
		std::string nrPositivesFilename = patOutPath + "/" + m_patientName + "_FastRipplesPerCh.txt";
		remove(nrPositivesFilename.c_str());
		std::ofstream nrPositivesFile;
		nrPositivesFile.open(nrPositivesFilename);
		nrPositivesFile << "nrMarks" << "\t" << "channelName" << "\n";
		long nrTotalPositives = 0;

		//#pragma omp parallel for reduction(+ : nrTotalPositives)
		for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
			int nrMarks = 0;
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double mark = epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx];
				if (mark > 0) {
					nrMarks++;
				}
			}
			if (nrMarks > 0) {
				chnnWithMarks[channIdx] = true;
				nrTotalPositives += nrMarks;
				std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
				channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
				nrPositivesFile << nrMarks << "\t" << channelName << "\n";
			}

		}
		nrPositivesFile << nrTotalPositives << "\t" << "TotalPositives" << "\n";
		nrPositivesFile.close();
	}

#pragma omp parallel for
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		if (chnnWithMarks[channIdx]) {
			std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
			channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
			int idNumber = m_useMontage ? m_montageLabels[channIdx].montageMOSSDET_Nr : m_unipolarLabels[channIdx].contactGlobalIdx;
			//std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization_Ripple.txt";
			std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization.txt";
			filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());

			std::ofstream trainingFile;
			if (!fileExists(filename.c_str())) {
				trainingFile.open(filename);
				//Header
				trainingFile << "Data \t" << "CH \t" << "Time \t";
				for (unsigned i = 0; i < epochFeatures.featureNames.size(); ++i) {
					trainingFile << epochFeatures.featureNames[i] << "\t";
				}
				trainingFile << "Output \n";
			}
			else {
				trainingFile.open(filename, std::ios::app);    // open file for appending
			}

			trainingFile.precision(32);
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			//Generate Excel readable File
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double sampleTime = epochFeatures.m_time[sampleIdx];
				unsigned currSample = epochFeatures.m_time[sampleIdx] * m_samplingRate;

				trainingFile << "d:\t" << channelName << "\t" << epochFeatures.m_time[sampleIdx] << "\t";

				unsigned nrFeats = epochFeatures.features.size();
				for (unsigned featIdx = 0; featIdx < nrFeats; ++featIdx) {
					unsigned featLength = epochFeatures.features[featIdx]->at(channIdx).size();
					if (featLength > 1 && sampleIdx < featLength) {
						trainingFile << epochFeatures.features[featIdx]->at(channIdx)[sampleIdx] << "\t";
					}
				}

				trainingFile << epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx] << "\n";
			}
			trainingFile.close();
		}
	}

	epochFeatures.normalizeFeatures();
	std::vector<int> featsToExclude;

	// bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes
	CorrelationFeatureSelection featureSelectorFastRipples(m_nrChannels, patOutPath, m_patientName);
	for (int i = 0; i <= 7; ++i) {
		featsToExclude.push_back(i);
	}
	featsToExclude.push_back(11);
	for (int i = 23; i <= 25; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 37; i <= 53; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 65; i <= 67; ++i) {
		featsToExclude.push_back(i);
	}
	featureSelectorFastRipples.setFeatsToExclude(featsToExclude);
	featureSelectorFastRipples.getEOI_MarkFeatureCorrelation(false, true, false, false, false, false, epochFeatures);	//Ripples

	return true;
}

bool HFO_Detector::generateEpochCharacterizationSpikes(EOI_Features &epochFeatures) {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patOutPath = outputPath + "/" + m_patientName;
	std::string patTrainingFilesPath = outputPath + "/" + m_patientName + "/TrainingFiles";

	_mkdir(outputPath.c_str());
	_mkdir(patOutPath.c_str());
	_mkdir(patTrainingFilesPath.c_str());


	std::vector<bool> chnnWithMarks(m_nrChannels);
	bool ripple = true;
	bool fastRipple = false;
	if (ripple) {
		getEpochMarkCoOccurrence(epochFeatures, SPIKE_LABEL);
		std::string nrPositivesFilename = patOutPath + "/" + m_patientName + "_SpikesPerCh.txt";
		remove(nrPositivesFilename.c_str());
		std::ofstream nrPositivesFile;
		nrPositivesFile.open(nrPositivesFilename);
		nrPositivesFile << "nrMarks" << "\t" << "channelName" << "\n";
		long nrTotalPositives = 0;

		//#pragma omp parallel for reduction(+ : nrTotalPositives)
		for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
			int nrMarks = 0;
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double mark = epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx];
				if (mark > 0) {
					nrMarks++;
				}
			}
			if (nrMarks > 0) {
				chnnWithMarks[channIdx] = true;
				nrTotalPositives += nrMarks;
				std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
				channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
				nrPositivesFile << nrMarks << "\t" << channelName << "\n";
			}

		}
		nrPositivesFile << nrTotalPositives << "\t" << "TotalPositives" << "\n";
		nrPositivesFile.close();
	}

#pragma omp parallel for
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		if (chnnWithMarks[channIdx]) {
			std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;
			channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
			int idNumber = m_useMontage ? m_montageLabels[channIdx].montageMOSSDET_Nr : m_unipolarLabels[channIdx].contactGlobalIdx;
			//std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization_Ripple.txt";
			std::string filename = patTrainingFilesPath + "/" + m_patientName + "_" + channelName + "_epochCharacterization.txt";
			filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());

			std::ofstream trainingFile;
			if (!fileExists(filename.c_str())) {
				trainingFile.open(filename);
				//Header
				trainingFile << "Data \t" << "CH \t" << "Time \t";
				for (unsigned i = 0; i < epochFeatures.featureNames.size(); ++i) {
					trainingFile << epochFeatures.featureNames[i] << "\t";
				}
				trainingFile << "Output \n";
			}
			else {
				trainingFile.open(filename, std::ios::app);    // open file for appending
			}

			trainingFile.precision(32);
			unsigned nrEpochs = epochFeatures.features[0]->at(0).size();
			//Generate Excel readable File
			for (unsigned sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double sampleTime = epochFeatures.m_time[sampleIdx];
				unsigned currSample = epochFeatures.m_time[sampleIdx] * m_samplingRate;

				trainingFile << "d:\t" << channelName << "\t" << epochFeatures.m_time[sampleIdx] << "\t";

				unsigned nrFeats = epochFeatures.features.size();
				for (unsigned featIdx = 0; featIdx < nrFeats; ++featIdx) {
					unsigned featLength = epochFeatures.features[featIdx]->at(channIdx).size();
					if (featLength > 1 && sampleIdx < featLength) {
						trainingFile << epochFeatures.features[featIdx]->at(channIdx)[sampleIdx] << "\t";
					}
				}

				trainingFile << epochFeatures.m_EOI_Mark.at(channIdx)[sampleIdx] << "\n";
			}
			trainingFile.close();
		}
	}

	return true;
	epochFeatures.normalizeFeatures();
	std::vector<int> featsToExclude;

	// bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes
	CorrelationFeatureSelection featureSelectorFastRipples(m_nrChannels, patOutPath, m_patientName);
	for (int i = 0; i <= 11; ++i) {
		featsToExclude.push_back(i);
	}
	for (int i = 40; i <= 67; ++i) {
		featsToExclude.push_back(i);
	}

	featureSelectorFastRipples.setFeatsToExclude(featsToExclude);
	featureSelectorFastRipples.getEOI_MarkFeatureCorrelation(false, false, false, false, false, true, epochFeatures);	//Ripples

	return true;
}

bool HFO_Detector::generateEpochCharacterizationSpikesHFO(EOI_Features &epochFeatures) {

	unsigned nrFeatures = epochFeatures.features.size();
	unsigned nrEpochs = epochFeatures.m_time.size();

	getEpochMarkCoOccurrence(epochFeatures, RIPPLE_LABEL);
	/***********************************************Copy Only Epochs with Spikes************************************************************************************/
	EOI_Features signalFeaturesSpikes;
	signalFeaturesSpikes.setMatrices(1);
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		for (long sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
			
			double epochLabel = epochFeatures.m_EOI_Mark[channIdx][sampleIdx];
			if (epochLabel >= SPIKE_LABEL - 0.1) {
				for (long featIdx = 0; featIdx < nrFeatures; ++featIdx) {
					double featVal = epochFeatures.features[featIdx]->at(channIdx)[sampleIdx];
					signalFeaturesSpikes.features.at(featIdx)->at(0).push_back(featVal);
				}
				signalFeaturesSpikes.m_EOI_Mark[0].push_back(epochLabel);
				double epochTime = epochFeatures.m_time[sampleIdx];
				signalFeaturesSpikes.m_time.push_back(epochTime);
			}
		}
	}

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patOutPath = outputPath + "/" + m_patientName;
	std::string patTrainingFilesPath = outputPath + "/" + m_patientName + "/TrainingFiles";


	_mkdir(outputPath.c_str());
	_mkdir(patOutPath.c_str());
	_mkdir(patTrainingFilesPath.c_str());

	std::string filename = patTrainingFilesPath + "/" + m_patientName + "_spike_vs_spikeHFO_EpochCharacterization.txt";
	filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());

	std::ofstream trainingFile;
	if (!fileExists(filename.c_str())) {
		trainingFile.open(filename);
		//Header
		trainingFile << "Data \t" << "CH \t" << "Time \t";
		for (unsigned i = 0; i < signalFeaturesSpikes.featureNames.size(); ++i) {
			trainingFile << signalFeaturesSpikes.featureNames[i] << "\t";
		}
		trainingFile << "Output \n";
	}
	else {
		trainingFile.open(filename, std::ios::app);    // open file for appending
	}

	trainingFile.precision(32);

	unsigned nrSpikeEpochs = signalFeaturesSpikes.features[0]->at(0).size();
	//Generate Excel readable File
	for (unsigned sampleIdx = 0; sampleIdx < nrSpikeEpochs; ++sampleIdx) {

		double sampleTime = signalFeaturesSpikes.m_time[sampleIdx];
		unsigned currSample = signalFeaturesSpikes.m_time[sampleIdx] * m_samplingRate;

		trainingFile << "d:\t" << "Spike_vs_SpikeHFO" << "\t" << signalFeaturesSpikes.m_time[sampleIdx] << "\t";

		unsigned nrFeats = signalFeaturesSpikes.features.size();
		for (unsigned featIdx = 0; featIdx < nrFeats; ++featIdx) {
			unsigned featLength = signalFeaturesSpikes.features[featIdx]->at(0).size();
			if (featLength > 1 && sampleIdx < featLength) {
				trainingFile << signalFeaturesSpikes.features[featIdx]->at(0)[sampleIdx] << "\t";
			}
		}

		trainingFile << signalFeaturesSpikes.m_EOI_Mark.at(0)[sampleIdx] << "\n";
	}
	trainingFile.close();



	signalFeaturesSpikes.normalizeFeatures();
	CorrelationFeatureSelection featureSelector(1, patOutPath, m_patientName);

	std::vector<int> featsToExclude;
	//Asymmetries
	featsToExclude.push_back(19);
	featsToExclude.push_back(33);
	featsToExclude.push_back(47);
	featsToExclude.push_back(61);
	//Spectral Peaks
	featsToExclude.push_back(2);
	featsToExclude.push_back(6);
	featsToExclude.push_back(10);

	// Oscills Ripple Range
	featsToExclude.push_back(44);
	featsToExclude.push_back(45);

	// Oscills FastRipple Range
	featsToExclude.push_back(58);
	featsToExclude.push_back(59);

	//FastRipple Features
	/*for (int i = 54; i <= 67; ++i) {
		featsToExclude.push_back(i);
	}*/
	featureSelector.setFeatsToExclude(featsToExclude);

	featureSelector.getEOI_MarkFeatureCorrelation(false, false, false, true, false, false, signalFeaturesSpikes);	//bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes

	//featureSelector.getEOI_MarkFeatureCorrelation(false, false, false, false, true, false, signalFeaturesSpikes);	//bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes

	return true;
}

void HFO_Detector::readChannelsFromFile(void) {
	
	m_useMontage ? m_selectedBipolarChanels.clear() : m_selectedUnipolarChanels.clear();

	//std::string fn = m_patientPath + m_patientName + "_AllChannels.txt";
	std::string fn = m_patientPath + m_patientName + "_EOI_ChannsFile.txt";
	//std::string fn = m_patientPath + m_patientName + "_RippleChannels.txt";
	//std::string fn = m_patientPath + m_patientName + "_FastRippleChannels.txt";
	//std::string fn = m_patientPath + m_patientName + "_SpikeChannels.txt";


	std::ifstream selectChanns;
	selectChanns.open(fn.c_str());

	std::string firstLine;
	std::getline(selectChanns, firstLine);

	while (!selectChanns.eof()) {
		std::string line;
		std::getline(selectChanns, line);
		std::stringstream ss(line);

		int channNr;
		std::string channName;
		ss >> channNr;
		ss >> channName;
		if (!channName.empty()) {
			m_useMontage ? m_selectedBipolarChanels.push_back(channName) : m_selectedUnipolarChanels.push_back(channName);
		}
	}

	selectChanns.close();
}
std::string HFO_Detector::ExePath() {
	char buff[FILENAME_MAX];
	_getcwd(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

void HFO_Detector::getSignalParams(double &samplingRate, double &analysisDuration, int &nrChannels) {
	samplingRate = m_samplingRate;
	analysisDuration = m_totalSampleCount / m_samplingRate;
	nrChannels = m_nrChannels;
}

void HFO_Detector::readAndWriteEOI_VisualAnnotations() {

	std::vector<EOI_Event> readAllAnnotations;
	long nrAnnotations = 0;

	/*****************************************************************************Read all annotations******************************************************************************************/
	std::string readFilename = m_patientPath + m_patientName + "_allAnnotations.txt";
	//std::string readFilename = m_patientPath + m_patientName + "_mixed_HFOAnnots_SpikeDetects_File.txt";

	std::ifstream eoiAnnotations;
	eoiAnnotations.open(readFilename.c_str());
	eoiAnnotations.precision(32);
	std::string firstLine;
	std::getline(eoiAnnotations, firstLine);
	while (!eoiAnnotations.eof()) {
		std::string line;
		std::getline(eoiAnnotations, line);
		std::stringstream ss(line);

		double sampleStartAndEnd;
		double duration;
		int refChannNr;
		std::string refChannName;

		EOI_Event eoiMark;
		/*ss >> eoiMark.channelName;
		ss >> eoiMark.description;
		ss >> refChannName;
		ss >> sampleStartAndEnd;
		ss >> duration;
		eoiMark.description += refChannName;
		eoiMark.description.erase(std::remove_if(eoiMark.description.begin(), eoiMark.description.end(), ::isspace), eoiMark.description.end());
		eoiMark.startTime = sampleStartAndEnd / m_samplingRate;
		eoiMark.endTime = sampleStartAndEnd / m_samplingRate + duration / m_samplingRate;*/

		/*ss >> eoiMark.channelName;
		ss >> eoiMark.channelNr;
		ss >> eoiMark.description;;
		ss >> refChannNr;
		ss >> sampleStartAndEnd;
		ss >> duration;
		ss >> eoiMark.startTime;
		ss >> eoiMark.endTime;*/

		ss >> eoiMark.description;
		ss >> eoiMark.channelName;
		ss >> eoiMark.channelNr;
		ss >> eoiMark.lowNeighborChannelName;
		ss >> eoiMark.lowNeighborChannelNr;
		ss >> eoiMark.highNeighborChannelName;
		ss >> eoiMark.highNeighborChannelNr;
		ss >> sampleStartAndEnd;
		ss >> sampleStartAndEnd;
		ss >> eoiMark.startTime;
		ss >> eoiMark.endTime;

		if (!eoiMark.channelName.empty() && eoiMark.startTime != 0 && eoiMark.endTime != 0) {
			if (getChannelNumber(eoiMark.channelName, eoiMark.channelNr)) {
				getChannelNumber(eoiMark.lowNeighborChannelName, eoiMark.lowNeighborChannelNr);
				getChannelNumber(eoiMark.highNeighborChannelName, eoiMark.highNeighborChannelNr);
				readAllAnnotations.push_back(eoiMark);
				nrAnnotations++;
			}
		}
	}
	
	eoiAnnotations.close();
	
	/*****************************************************************************Write all annotations******************************************************************************************/
	m_allMarks.resize(m_nrChannels);
	m_markedEOI.resize(m_nrChannels);
	m_markedSpikesHFO.resize(m_nrChannels);
	m_markedSpikesAlone.resize(m_nrChannels);
	m_markedAllSpikes.resize(m_nrChannels);

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName + "/";
	_mkdir(outputPath.c_str());
	_mkdir(patPath.c_str());
	std::string writeFilename = patPath + m_patientName + "_allAnnotations.txt";
	remove(writeFilename.c_str());

	std::ofstream eoiAnnotWriteFile;
	if (!fileExists(writeFilename.c_str())) {
		eoiAnnotWriteFile.open(writeFilename);
		//Header
		eoiAnnotWriteFile << "Mark\t" << "Channel\t" << "ChannelNr\t" << "LowerNeighborChName\t" << "LowerNeighborChNr\t" << "UpperNeighborChName\t" << "UpperNeighborChNr\t" << "StartSample\t" << "EndSample\t" << "StartTime\t" << "EndTime\n";
	}

	for (int annotIdx = 0; annotIdx < nrAnnotations; ++annotIdx) {

		std::string label = readAllAnnotations[annotIdx].description;
		bool foundFastRipple = label.find("ast") != std::string::npos;
		bool foundRipple = ((label.find("ipple") != std::string::npos) || (label.find("HFO") != std::string::npos)) && !foundFastRipple;
		bool foundSpike = label.find("ike") != std::string::npos;


		long annotChannNr = readAllAnnotations[annotIdx].channelNr;
		long lowerNeighborChNr = annotChannNr;
		long upperNeighborChNr = annotChannNr; 

		
		if (m_useMontage) {
			getNeighboringNumber(annotChannNr, lowerNeighborChNr, upperNeighborChNr);
			readAllAnnotations[annotIdx].lowNeighborChannelNr = lowerNeighborChNr;
			readAllAnnotations[annotIdx].highNeighborChannelNr = upperNeighborChNr;

			if (readAllAnnotations[annotIdx].lowNeighborChannelNr > readAllAnnotations[annotIdx].highNeighborChannelNr)
				while(1);
		}

		m_allMarks[annotChannNr].push_back(readAllAnnotations[annotIdx]);
		if (foundRipple || foundFastRipple) {
			
			double duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;
			double center = readAllAnnotations[annotIdx].startTime + duration / 2;
			//readAllAnnotations[annotIdx].startTime = center - 0.025;
			//readAllAnnotations[annotIdx].endTime = center + 0.025;
			readAllAnnotations[annotIdx].duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;

			m_markedEOI[annotChannNr].push_back(readAllAnnotations[annotIdx]);

			if (readAllAnnotations[annotIdx].startTime < m_firstMarkedSample)
				m_firstMarkedSample = readAllAnnotations[annotIdx].startTime;
			if (readAllAnnotations[annotIdx].endTime > m_lastMarkedSample)
				m_lastMarkedSample = readAllAnnotations[annotIdx].endTime;
		}
		else if(foundSpike) {
			double duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;
			double center = readAllAnnotations[annotIdx].startTime + duration / 2;
			//readAllAnnotations[annotIdx].startTime = center - 0.1;
			//readAllAnnotations[annotIdx].endTime = center + 0.1;
			readAllAnnotations[annotIdx].duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;
			m_markedSpikesAlone[annotChannNr].push_back(readAllAnnotations[annotIdx]);
		}

		if (foundSpike) {
			double duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;
			double center = readAllAnnotations[annotIdx].startTime + duration / 2;
			//readAllAnnotations[annotIdx].startTime = center - 0.1;
			//readAllAnnotations[annotIdx].endTime = center + 0.1;
			readAllAnnotations[annotIdx].duration = readAllAnnotations[annotIdx].endTime - readAllAnnotations[annotIdx].startTime;
			m_markedAllSpikes[annotChannNr].push_back(readAllAnnotations[annotIdx]);
		}

		eoiAnnotWriteFile << readAllAnnotations[annotIdx].description << "\t"
							<< readAllAnnotations[annotIdx].channelName << "\t" << readAllAnnotations[annotIdx].channelNr << "\t";											//Mark Channel
		if (m_useMontage) {
			eoiAnnotWriteFile << m_montageLabels[lowerNeighborChNr].montageName << "\t" << m_montageLabels[lowerNeighborChNr].montageMOSSDET_Nr << "\t"						//Lower Neighbor Channel
							<< m_montageLabels[upperNeighborChNr].montageName << "\t" << m_montageLabels[upperNeighborChNr].montageMOSSDET_Nr << "\t";						//Upper Neighbor Channel
		}
		else {
			eoiAnnotWriteFile << readAllAnnotations[annotIdx].channelName << "\t" << readAllAnnotations[annotIdx].channelNr << "\t"											//Mark Channel
							<< readAllAnnotations[annotIdx].channelName << "\t" << readAllAnnotations[annotIdx].channelNr << "\t";											//Mark Channel
		}
		
		eoiAnnotWriteFile << (long) readAllAnnotations[annotIdx].startTime * m_samplingRate << "\t" << (long) readAllAnnotations[annotIdx].endTime * m_samplingRate << "\t"				// Start and end samples
						<< readAllAnnotations[annotIdx].startTime << "\t" << readAllAnnotations[annotIdx].endTime << "\n";												// Start and end times

	}
	eoiAnnotWriteFile.close();

	getSpikesAndSpikeHFOs();

	/*****************************************************************************Write HFO annotations******************************************************************************************/
	writeFilename = patPath + m_patientName + "_eoiAnnotations.txt";
	remove(writeFilename.c_str());

	std::ofstream allHFOAnnotWriteFile;
	if (!fileExists(writeFilename.c_str())) {
		allHFOAnnotWriteFile.open(writeFilename);
		//Header
		allHFOAnnotWriteFile << "Mark\t" << "Channel\t" << "ChannelNr\t" << "LowerNeighborChName\t" << "LowerNeighborChNr\t" << "UpperNeighborChName\t" << "UpperNeighborChNr\t" << "StartSample\t" << "EndSample\t" << "StartTime\t" << "EndTime\n";
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		
		long nrEOI_Annotations = m_markedEOI[channIdx].size();
		
		for (int annotIdx = 0; annotIdx < nrEOI_Annotations; ++annotIdx) {
			long annotChannNr = m_markedEOI[channIdx][annotIdx].channelNr;
			long upperNeighborChNr = annotChannNr;
			long lowerNeighborChNr = annotChannNr;

			if (m_useMontage) {
				getNeighboringNumber(annotChannNr, lowerNeighborChNr, upperNeighborChNr);
				if (lowerNeighborChNr > upperNeighborChNr)
					while (1);
			}

			allHFOAnnotWriteFile << m_markedEOI[channIdx][annotIdx].description << "\t"
				<< m_markedEOI[channIdx][annotIdx].channelName << "\t" << m_markedEOI[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			if (m_useMontage) {
				allHFOAnnotWriteFile << m_montageLabels[lowerNeighborChNr].montageName << "\t" << m_montageLabels[lowerNeighborChNr].montageMOSSDET_Nr << "\t"						//Lower Neighbor Channel
					<< m_montageLabels[upperNeighborChNr].montageName << "\t" << m_montageLabels[upperNeighborChNr].montageMOSSDET_Nr << "\t";						//Upper Neighbor Channel
			}
			else {
				allHFOAnnotWriteFile << m_markedEOI[channIdx][annotIdx].channelName << "\t" << m_markedEOI[channIdx][annotIdx].channelNr << "\t"											//Mark Channel
					<< m_markedEOI[channIdx][annotIdx].channelName << "\t" << m_markedEOI[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			}

			allHFOAnnotWriteFile << (long) m_markedEOI[channIdx][annotIdx].startTime * m_samplingRate << "\t" << (long) m_markedEOI[channIdx][annotIdx].endTime * m_samplingRate << "\t"				// Start and end samples
				<< m_markedEOI[channIdx][annotIdx].startTime << "\t" << m_markedEOI[channIdx][annotIdx].endTime << "\n";												// Start and end times
		}
	}
	allHFOAnnotWriteFile.close();

	/*****************************************************************************Write Spike+HFO annotations******************************************************************************************/
	writeFilename = patPath + m_patientName + "_spikeHFOAnnotations.txt";
	remove(writeFilename.c_str());

	std::ofstream allSpikeHFOAnnotWriteFile;
	if (!fileExists(writeFilename.c_str())) {
		allSpikeHFOAnnotWriteFile.open(writeFilename);
		//Header
		allSpikeHFOAnnotWriteFile << "Mark\t" << "Channel\t" << "ChannelNr\t" << "LowerNeighborChName\t" << "LowerNeighborChNr\t" << "UpperNeighborChName\t" << "UpperNeighborChNr\t" << "StartSample\t" << "EndSample\t" << "StartTime\t" << "EndTime\n";
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {

		long nrEOI_Annotations = m_markedSpikesHFO[channIdx].size();

		for (int annotIdx = 0; annotIdx < nrEOI_Annotations; ++annotIdx) {
			long annotChannNr = m_markedSpikesHFO[channIdx][annotIdx].channelNr;
			long upperNeighborChNr = annotChannNr;
			long lowerNeighborChNr = annotChannNr;

			if (m_useMontage) {
				getNeighboringNumber(annotChannNr, lowerNeighborChNr, upperNeighborChNr);
				if (lowerNeighborChNr > upperNeighborChNr)
					while (1);
			}

			allSpikeHFOAnnotWriteFile << m_markedSpikesHFO[channIdx][annotIdx].description << "\t"
				<< m_markedSpikesHFO[channIdx][annotIdx].channelName << "\t" << m_markedSpikesHFO[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			if (m_useMontage) {
				allSpikeHFOAnnotWriteFile << m_montageLabels[lowerNeighborChNr].montageName << "\t" << m_montageLabels[lowerNeighborChNr].montageMOSSDET_Nr << "\t"						//Lower Neighbor Channel
					<< m_montageLabels[upperNeighborChNr].montageName << "\t" << m_montageLabels[upperNeighborChNr].montageMOSSDET_Nr << "\t";						//Upper Neighbor Channel
			}
			else {
				allSpikeHFOAnnotWriteFile << m_markedSpikesHFO[channIdx][annotIdx].channelName << "\t" << m_markedSpikesHFO[channIdx][annotIdx].channelNr << "\t"											//Mark Channel
					<< m_markedSpikesHFO[channIdx][annotIdx].channelName << "\t" << m_markedSpikesHFO[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			}

			allSpikeHFOAnnotWriteFile << (long) m_markedSpikesHFO[channIdx][annotIdx].startTime * m_samplingRate << "\t" << (long) m_markedSpikesHFO[channIdx][annotIdx].endTime * m_samplingRate << "\t"				// Start and end samples
				<< m_markedSpikesHFO[channIdx][annotIdx].startTime << "\t" << m_markedSpikesHFO[channIdx][annotIdx].endTime << "\n";												// Start and end times
		}
	}
	allSpikeHFOAnnotWriteFile.close();

	/*****************************************************************************Write only Spike annotations******************************************************************************************/
	writeFilename = patPath + m_patientName + "_spikeOnlyAnnotations.txt";
	remove(writeFilename.c_str());

	std::ofstream spikeOnlyAnnotWriteFile;
	if (!fileExists(writeFilename.c_str())) {
		spikeOnlyAnnotWriteFile.open(writeFilename);
		//Header
		spikeOnlyAnnotWriteFile << "Mark\t" << "Channel\t" << "ChannelNr\t" << "LowerNeighborChName\t" << "LowerNeighborChNr\t" << "UpperNeighborChName\t" << "UpperNeighborChNr\t" << "StartSample\t" << "EndSample\t" << "StartTime\t" << "EndTime\n";
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {

		long nrEOI_Annotations = m_markedSpikesAlone[channIdx].size();

		for (int annotIdx = 0; annotIdx < nrEOI_Annotations; ++annotIdx) {
			long annotChannNr = m_markedSpikesAlone[channIdx][annotIdx].channelNr;
			long upperNeighborChNr = annotChannNr;
			long lowerNeighborChNr = annotChannNr;

			if (m_useMontage) {
				getNeighboringNumber(annotChannNr, lowerNeighborChNr, upperNeighborChNr);
				if (lowerNeighborChNr > upperNeighborChNr)
					while (1);
			}

			spikeOnlyAnnotWriteFile << m_markedSpikesAlone[channIdx][annotIdx].description << "\t"
				<< m_markedSpikesAlone[channIdx][annotIdx].channelName << "\t" << m_markedSpikesAlone[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			if (m_useMontage) {
				spikeOnlyAnnotWriteFile << m_montageLabels[lowerNeighborChNr].montageName << "\t" << m_montageLabels[lowerNeighborChNr].montageMOSSDET_Nr << "\t"						//Lower Neighbor Channel
					<< m_montageLabels[upperNeighborChNr].montageName << "\t" << m_montageLabels[upperNeighborChNr].montageMOSSDET_Nr << "\t";						//Upper Neighbor Channel
			}
			else {
				spikeOnlyAnnotWriteFile << m_markedSpikesAlone[channIdx][annotIdx].channelName << "\t" << m_markedSpikesAlone[channIdx][annotIdx].channelNr << "\t"											//Mark Channel
					<< m_markedSpikesAlone[channIdx][annotIdx].channelName << "\t" << m_markedSpikesAlone[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			}

			spikeOnlyAnnotWriteFile << (long)m_markedSpikesAlone[channIdx][annotIdx].startTime * m_samplingRate << "\t" << (long)m_markedSpikesAlone[channIdx][annotIdx].endTime * m_samplingRate << "\t"				// Start and end samples
				<< m_markedSpikesAlone[channIdx][annotIdx].startTime << "\t" << m_markedSpikesAlone[channIdx][annotIdx].endTime << "\n";												// Start and end times
		}
	}
	spikeOnlyAnnotWriteFile.close();

	/*****************************************************************************Write all Spike annotations******************************************************************************************/
	writeFilename = patPath + m_patientName + "_allSpikeAnnotations.txt";
	remove(writeFilename.c_str());

	std::ofstream allSpikeAnnotWriteFile;
	if (!fileExists(writeFilename.c_str())) {
		allSpikeAnnotWriteFile.open(writeFilename);
		//Header
		allSpikeAnnotWriteFile << "Mark\t" << "Channel\t" << "ChannelNr\t" << "LowerNeighborChName\t" << "LowerNeighborChNr\t" << "UpperNeighborChName\t" << "UpperNeighborChNr\t" << "StartSample\t" << "EndSample\t" << "StartTime\t" << "EndTime\n";
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {

		long nrEOI_Annotations = m_markedAllSpikes[channIdx].size();

		for (int annotIdx = 0; annotIdx < nrEOI_Annotations; ++annotIdx) {
			long annotChannNr = m_markedAllSpikes[channIdx][annotIdx].channelNr;
			long upperNeighborChNr = annotChannNr;
			long lowerNeighborChNr = annotChannNr;

			if (m_useMontage) {
				getNeighboringNumber(annotChannNr, lowerNeighborChNr, upperNeighborChNr);
				if (lowerNeighborChNr > upperNeighborChNr)
					while (1);
			}

			allSpikeAnnotWriteFile << m_markedAllSpikes[channIdx][annotIdx].description << "\t"
				<< m_markedAllSpikes[channIdx][annotIdx].channelName << "\t" << m_markedAllSpikes[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			if (m_useMontage) {
				allSpikeAnnotWriteFile << m_montageLabels[lowerNeighborChNr].montageName << "\t" << m_montageLabels[lowerNeighborChNr].montageMOSSDET_Nr << "\t"						//Lower Neighbor Channel
					<< m_montageLabels[upperNeighborChNr].montageName << "\t" << m_montageLabels[upperNeighborChNr].montageMOSSDET_Nr << "\t";						//Upper Neighbor Channel
			}
			else {
				allSpikeAnnotWriteFile << m_markedAllSpikes[channIdx][annotIdx].channelName << "\t" << m_markedAllSpikes[channIdx][annotIdx].channelNr << "\t"											//Mark Channel
					<< m_markedAllSpikes[channIdx][annotIdx].channelName << "\t" << m_markedAllSpikes[channIdx][annotIdx].channelNr << "\t";											//Mark Channel
			}

			allSpikeAnnotWriteFile << (long)m_markedAllSpikes[channIdx][annotIdx].startTime * m_samplingRate << "\t" << (long)m_markedAllSpikes[channIdx][annotIdx].endTime * m_samplingRate << "\t"				// Start and end samples
				<< m_markedAllSpikes[channIdx][annotIdx].startTime << "\t" << m_markedAllSpikes[channIdx][annotIdx].endTime << "\n";												// Start and end times
		}
	}
	allSpikeAnnotWriteFile.close();
}

bool HFO_Detector::getChannelNumber(std::string markChannelName, unsigned int &channNr) {

//#pragma omp parallel for
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {

		std::string channelName = m_useMontage ? m_montageLabels[channIdx].montageName : m_unipolarLabels[channIdx].contactName;

		markChannelName.erase(std::remove_if(markChannelName.begin(), markChannelName.end(), ::isspace), markChannelName.end());
		channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());

		markChannelName.erase(std::remove(markChannelName.begin(), markChannelName.end(), '-'), markChannelName.end());
		channelName.erase(std::remove(channelName.begin(), channelName.end(), '-'), channelName.end());

		markChannelName.erase(std::remove(markChannelName.begin(), markChannelName.end(), '_'), markChannelName.end());
		channelName.erase(std::remove(channelName.begin(), channelName.end(), '_'), channelName.end());

		markChannelName.erase(std::remove(markChannelName.begin(), markChannelName.end(), '#'), markChannelName.end());
		channelName.erase(std::remove(channelName.begin(), channelName.end(), '#'), channelName.end());

		int compareResult = markChannelName.compare(channelName);
		if (compareResult == 0) {
			channNr = channIdx;
			return true;
		}
		
	}

	return false;
}

/*bool HFO_Detector::getNeighboringNumber(long annotChannNr, long &lowerNeighborChNr, long &upperNeighborChNr) {

	bool sameElectrodeUp, sameElectrodeDown = false;
	if (annotChannNr > 0) {
		sameElectrodeDown = (m_montageLabels[annotChannNr].secondElectrodeName.compare(m_montageLabels[annotChannNr - 1].firstElectrodeName) == 0) && (m_montageLabels[annotChannNr].firstElectrodeName.compare(m_montageLabels[annotChannNr - 1].secondElectrodeName) == 0);
		if (sameElectrodeDown) {
			lowerNeighborChNr = annotChannNr - 1;
		}

		if (annotChannNr < m_montageLabels.size() - 1) {
			sameElectrodeUp = (m_montageLabels[annotChannNr].secondElectrodeName.compare(m_montageLabels[annotChannNr + 1].firstElectrodeName) == 0) && (m_montageLabels[annotChannNr].firstElectrodeName.compare(m_montageLabels[annotChannNr + 1].secondElectrodeName) == 0);
			if (sameElectrodeUp) {
				upperNeighborChNr = annotChannNr + 1;
			}
		}
		if (lowerNeighborChNr > upperNeighborChNr)
			while (1);
		return sameElectrodeUp && sameElectrodeDown;

	}
	else if(annotChannNr < m_montageLabels.size()-1){
		sameElectrodeUp = (m_montageLabels[annotChannNr].secondElectrodeName.compare(m_montageLabels[annotChannNr + 1].firstElectrodeName) == 0) && (m_montageLabels[annotChannNr].firstElectrodeName.compare(m_montageLabels[annotChannNr + 1].secondElectrodeName) == 0);
		if (sameElectrodeUp) {
			upperNeighborChNr = annotChannNr + 1;
			if (lowerNeighborChNr > upperNeighborChNr)
				while (1);
			return true;
		}
	}

	return false;
}*/

bool HFO_Detector::getNeighboringNumber(long annotChannNr, long &lowerNeighborChNr, long &upperNeighborChNr) {

	//getLowerNeighbor
	std::string analyzedChann = m_montageLabels[annotChannNr].montageName;
	std::string possibleNeighborName = m_montageLabels[annotChannNr].firstElectrodeName + std::to_string(m_montageLabels[annotChannNr].firstContactNr - 1) + m_montageLabels[annotChannNr].secondElectrodeName + std::to_string(m_montageLabels[annotChannNr].firstContactNr);
	//std::string possibleNeighborName = m_montageLabels[annotChannNr].firstElectrodeName + std::to_string(m_montageLabels[annotChannNr].firstContactNr) + m_montageLabels[annotChannNr].secondElectrodeName + std::to_string(m_montageLabels[annotChannNr].firstContactNr + 1);
	
	bool neighboringChannelFound = false;
	for (int i = 0; i < m_montageLabels.size(); ++i) {
		std::string currentChannel = m_montageLabels[i].montageName;
		currentChannel.erase(std::remove_if(currentChannel.begin(), currentChannel.end(), ::isspace), currentChannel.end());
		possibleNeighborName.erase(std::remove_if(possibleNeighborName.begin(), possibleNeighborName.end(), ::isspace), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '-'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '-'), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '_'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '_'), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '#'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '#'), possibleNeighborName.end());

		int compareResult = currentChannel.compare(possibleNeighborName);
		if (compareResult == 0) {
			lowerNeighborChNr = i;
			neighboringChannelFound = true;
			break;
		}
	}
	if (!neighboringChannelFound) {
		lowerNeighborChNr = annotChannNr;
	}

	//getHigherNeighbor
	possibleNeighborName = m_montageLabels[annotChannNr].firstElectrodeName + std::to_string(m_montageLabels[annotChannNr].secondContactNr) + m_montageLabels[annotChannNr].secondElectrodeName + std::to_string(m_montageLabels[annotChannNr].secondContactNr + 1);
	neighboringChannelFound = false;
	for (int i = 0; i < m_montageLabels.size(); ++i) {
		std::string currentChannel = m_montageLabels[i].montageName;
		currentChannel.erase(std::remove_if(currentChannel.begin(), currentChannel.end(), ::isspace), currentChannel.end());
		possibleNeighborName.erase(std::remove_if(possibleNeighborName.begin(), possibleNeighborName.end(), ::isspace), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '-'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '-'), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '_'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '_'), possibleNeighborName.end());

		currentChannel.erase(std::remove(currentChannel.begin(), currentChannel.end(), '#'), currentChannel.end());
		possibleNeighborName.erase(std::remove(possibleNeighborName.begin(), possibleNeighborName.end(), '#'), possibleNeighborName.end());

		int compareResult = currentChannel.compare(possibleNeighborName);
		if (compareResult == 0) {
			upperNeighborChNr = i;
			neighboringChannelFound = true;
			break;
		}
	}
	if (!neighboringChannelFound) {
		upperNeighborChNr = annotChannNr;
	}
	
	return false;
}

bool HFO_Detector::getEpochMarkCoOccurrence(EOI_Features &epochFeatures, int eventType) {
	
	double epochTimeLength = ((double) ceil(EPOCH * m_samplingRate)) / m_samplingRate;
	double minSharedPercent = 0.80;
	
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		unsigned nrEpochs = epochFeatures.m_time.size();
#pragma omp parallel for
		for (long sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
			epochFeatures.m_EOI_Mark[channIdx][sampleIdx] = 0;
		}
	}

	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
		unsigned nrEpochs = epochFeatures.m_time.size();

		for (long sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
			double epochStart = epochFeatures.m_time[sampleIdx]- epochTimeLength;
			double epochEnd = epochFeatures.m_time[sampleIdx];
			double epochLabel = 0.0;
			epochFeatures.m_EOI_Mark[channIdx][sampleIdx] = 0;

			long nrAnnotations = m_markedEOI[channIdx].size();
			for (int annotIdx = 0; annotIdx < nrAnnotations; ++annotIdx) {
				std::string markDescription = m_markedEOI[channIdx][annotIdx].description;
				double markStart = m_markedEOI[channIdx][annotIdx].startTime;
				double markEnd = m_markedEOI[channIdx][annotIdx].endTime;
				long markLowNeighborChannNr = m_markedEOI[channIdx][annotIdx].lowNeighborChannelNr;
				long markHighNeighborChannNr = m_markedEOI[channIdx][annotIdx].highNeighborChannelNr;
				if (areEventsOverlapping(epochStart, epochEnd, markStart, markEnd, minSharedPercent)) {
					//if (eventType == getHFO_Code(m_markedEOI[channIdx][annotIdx].description)) {
						epochLabel = eventType;
					epochLabel = getHFO_Code(m_markedEOI[channIdx][annotIdx].description);

						if (CONSIDER_PROPAGATION_HFO) {						// label neighbor only of it hasn't been labelled yet
							if (channIdx != markLowNeighborChannNr) {
								if (compareDoublesSmaller(epochFeatures.m_EOI_Mark[markLowNeighborChannNr][sampleIdx], 0.1)) {
									epochFeatures.m_EOI_Mark[markLowNeighborChannNr][sampleIdx] = epochLabel + PROPAGATED_ANNOT_LABEL;
								}
							}
							if (channIdx != markHighNeighborChannNr) {
								if (compareDoublesSmaller(epochFeatures.m_EOI_Mark[markHighNeighborChannNr][sampleIdx], 0.1)) {
									epochFeatures.m_EOI_Mark[markHighNeighborChannNr][sampleIdx] = epochLabel + PROPAGATED_ANNOT_LABEL;
								}
							}
						}

						break;
					//}
				}
			}

			nrAnnotations = m_markedAllSpikes[channIdx].size();
			for (int annotIdx = 0; annotIdx < nrAnnotations; ++annotIdx) {
				double markStart = m_markedAllSpikes[channIdx][annotIdx].startTime;
				double markEnd = m_markedAllSpikes[channIdx][annotIdx].endTime;
				long markLowNeighborChannNr = m_markedAllSpikes[channIdx][annotIdx].lowNeighborChannelNr;
				long markHighNeighborChannNr = m_markedAllSpikes[channIdx][annotIdx].highNeighborChannelNr;
				if (areEventsOverlapping(epochStart, epochEnd, markStart, markEnd, minSharedPercent)) {
					epochLabel += SPIKE_LABEL;

					if (CONSIDER_PROPAGATION_SPIKES) {							// with spikes, all neighbors can be labelled as having a spike
						if (channIdx != markLowNeighborChannNr) {
							if (compareDoublesSmaller(epochFeatures.m_EOI_Mark[markLowNeighborChannNr][sampleIdx], SPIKE_LABEL))
								epochFeatures.m_EOI_Mark[markLowNeighborChannNr][sampleIdx] += SPIKE_LABEL;
						}
						if (channIdx != markHighNeighborChannNr) {
							if (compareDoublesSmaller(epochFeatures.m_EOI_Mark[markHighNeighborChannNr][sampleIdx], 0.0))
								epochFeatures.m_EOI_Mark[markHighNeighborChannNr][sampleIdx] += SPIKE_LABEL;
						}
					}

					break;
				}
			}

			epochFeatures.m_EOI_Mark[channIdx][sampleIdx] = epochLabel;
		}
	}
	return true;
}


bool HFO_Detector::areEventsOverlapping(double startTimeA, double endTimeA, double startTimeB, double endTimeB, double minSharedPercent) {
	bool overlap = false;

	double durationA = endTimeA - startTimeA;
	double durationB = endTimeB - startTimeB;

	double shared = 0;
	double minShared = durationA * minSharedPercent < durationB * minSharedPercent ? durationA * minSharedPercent : durationB * minSharedPercent;

	if (startTimeA >= startTimeB && endTimeA <= endTimeB) {			//Epoch completely within Mark
		shared = durationA;
		shared = 1;
	}
	else if (startTimeA <= startTimeB && endTimeA >= endTimeB) {		//Mark completely within Epoch
		shared = durationB;
		shared = 1;
	}
	else if (startTimeA >= startTimeB && startTimeA <= endTimeB) {		//Epoch beginning is within the Mark
		if (endTimeA >= endTimeB) {
			shared = endTimeB - startTimeA;
		}
		else {
			shared = durationA;
			shared = 1;
		}
	}
	else if (startTimeA <= startTimeB && endTimeA >= startTimeB) {		//Epoch end is within the Mark
		if (endTimeA <= endTimeB) {
			shared = endTimeA - startTimeB;
		}
		else {
			shared = durationB;
			shared = 1;
		}
	}

	return (shared >= minShared);
}
int HFO_Detector::getHFO_Code(std::string mark) {

	int eventCode = 0;

	bool foundFastRipple = mark.find("ast") != std::string::npos;
	bool foundRipple = ((mark.find("ipple") != std::string::npos) || (mark.find("HFO") != std::string::npos)) && !foundFastRipple;
	bool foundRippleAndFastRipple = false;
	if (mark.find("ipple") != std::string::npos) {
		mark.erase(mark.find("ipple"), 5);
		foundRippleAndFastRipple = mark.find("ipple") != std::string::npos;
		if (foundRippleAndFastRipple) {
			foundRipple = false;
			foundFastRipple = false;
		}
	}
	
	if (foundRipple) {
		eventCode = RIPPLE_LABEL;
	}
	else if (foundFastRipple) {
		eventCode = FAST_RIPPLE_LABEL;
	}
	else if (foundRippleAndFastRipple) {
		eventCode = RIPPLE_FR_LABEL;
	}

	return eventCode;
}

void HFO_Detector::getSpecificStartAndEndTimes() {

	if (m_patientName.find("0_NM") != std::string::npos) {
		m_analysisStartTime = 24.00;
		m_analysisEndTime = 84;
	}
	else if (m_patientName.find("0_rksegment1_Segment_1") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 61;
	}
	else if (m_patientName.find("0_rksegment2_Segment_1") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 179;
	}
	else if (m_patientName.find("09_KA") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 60;
	}
	else if (m_patientName.find("15_SchM_2000chan1_128_Segment_1") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 62;
	}
	else if (m_patientName.find("17_StV") != std::string::npos) {
		m_analysisStartTime = 125;
		m_analysisEndTime = 187;
	}
	else if (m_patientName.find("20_WoT") != std::string::npos) {
		m_analysisStartTime = 68;
		m_analysisEndTime = 130;
	}
	else if ((m_patientName.find("_SCHLAF_") != std::string::npos) || (m_patientName.find("_SCHLAF1_") != std::string::npos)) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 181;
	}
	else if (m_patientName.find("KUSSA") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 300;
	}
	else if (m_patientName.find("KONV") != std::string::npos || m_patientName.find("FERTIGMARK") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 62;
	} 
	else {
		double samplingRate;
		long firstAnalysisSample = 0, lastAnalysisSample = 0;
		std::string line, varName;

		std::string fn = m_patientPath + m_patientName + "_AnalysisTimeFile.txt";
		std::ifstream analysisTimeFile;
		analysisTimeFile.open(fn.c_str());

		std::getline(analysisTimeFile, line);
		std::stringstream first(line);
		first >> varName;
		first >> samplingRate;

		std::getline(analysisTimeFile, line);
		std::stringstream second(line);
		second >> varName;
		second >> firstAnalysisSample;
	
		std::getline(analysisTimeFile, line);
		std::stringstream third(line);
		third >> varName;
		third >> lastAnalysisSample;

		m_analysisStartTime = firstAnalysisSample / samplingRate;
		m_analysisEndTime = lastAnalysisSample / samplingRate;

		analysisTimeFile.close();
	}

	/*if (m_patientName.find("-SNR") != std::string::npos) {
		m_analysisStartTime = 0;
		m_analysisEndTime = 120;
	}*/
}

void HFO_Detector::getSpikesAndSpikeHFOs(void) {

	double minSharedPercent = 0.25;
		
//#pragma omp parallel for
	for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
			
		long nrSpike_Annotations = m_markedAllSpikes[channIdx].size();
			
		for (int spikeAnnotIdx = 0; spikeAnnotIdx < nrSpike_Annotations; ++spikeAnnotIdx) {
				
			double spikeStart = m_markedAllSpikes[channIdx][spikeAnnotIdx].startTime;
			double spikeEnd = m_markedAllSpikes[channIdx][spikeAnnotIdx].endTime;
			double  spikeDuration = spikeEnd - spikeStart;

			bool coincident = false;

			long nrHFO_Annotations = m_markedEOI[channIdx].size();

			for (int hfoAnnotIdx = 0; hfoAnnotIdx < nrHFO_Annotations; ++hfoAnnotIdx) {
					
				double hfoStart = m_markedEOI[channIdx][hfoAnnotIdx].startTime;
				double hfoEnd = m_markedEOI[channIdx][hfoAnnotIdx].endTime;
				double hfoDuration = hfoEnd - hfoStart;

				double shared = 0;
				double minShared = spikeDuration * minSharedPercent < hfoDuration * minSharedPercent ? spikeDuration * minSharedPercent : hfoDuration * minSharedPercent;

				if (spikeStart >= hfoStart && spikeEnd <= hfoEnd) {			//Epoch completely within Mark
					shared = spikeDuration;
				}
				else if (spikeStart <= hfoStart && spikeEnd >= hfoEnd) {		//Mark completely within Epoch
					shared = hfoDuration;
				}
				else if (spikeStart >= hfoStart && spikeStart <= hfoEnd) {		//Epoch beginning is within the Mark
					if (spikeEnd >= hfoEnd) {
						shared = hfoEnd - spikeStart;
					}
					else {
						shared = spikeDuration;
					}
				}
				else if (spikeStart <= hfoStart && spikeEnd >= hfoStart) {		//Epoch end is within the Mark
					if (spikeEnd <= hfoEnd) {
						shared = spikeEnd - hfoStart;
					}
					else {
						shared = hfoDuration;
					}
				}

				if (shared >= minShared) {
					coincident = true;
					m_markedEOI[channIdx][hfoAnnotIdx].description = "Spike+" + m_markedEOI[channIdx][hfoAnnotIdx].description;
					//break;
				}
			}

			if (coincident) {
				m_markedSpikesHFO[channIdx].push_back(m_markedAllSpikes[channIdx][spikeAnnotIdx]);
			}
			else {
				m_markedSpikesAlone[channIdx].push_back(m_markedAllSpikes[channIdx][spikeAnnotIdx]);
			}
		}
	}
}

				

void HFO_Detector::writeVectorToFile(int channNr, double time, std::vector<double> frequency, std::vector<double> vectorB) {
	
	std::string channelName = m_useMontage ? m_montageLabels[channNr].montageName : m_unipolarLabels[channNr].contactName;

	bool hippocampalElectrode = channNr == 0;// channelName.find("HC") != std::string::npos;
	if (hippocampalElectrode) {
		std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
		std::string patPath = outputPath + "/" + m_patientName + "/";
		std::string mwtPath = patPath + "/" + "MWT/";
		_mkdir(mwtPath.c_str());

		std::string filename = mwtPath + channelName +"SampleMWT_Plot_" + std::to_string(m_waveletParams->waveletOscillations) + "_Oscills.txt";

		std::ofstream vectorsFile;
		std::ifstream infile(filename);
		if (!infile.good()) {
			//Header
			vectorsFile.open(filename);
			vectorsFile << "Time" << "\t" << "Frequency" << "\t" << "Power" << "\n";
		}
		else {
			vectorsFile.open(filename, std::ios::app);    // open file for appending
		}

		for (long i = 0; i < vectorB.size(); ++i) {
			vectorsFile << time << "\t" << frequency[i] << "\t" << vectorB[i] << "\n";
		}

		vectorsFile.close();
	}
}

void HFO_Detector::readChannelsAndResection() {

	std::string fn = m_patientPath + m_patientName + "_channelsAndResectionFile.txt";

	if (fileExists(fn.c_str())) {

		std::ifstream channelsAndResectionFile;
		channelsAndResectionFile.open(fn.c_str());
		std::string firstLine;
		std::getline(channelsAndResectionFile, firstLine);

		while (!channelsAndResectionFile.eof()) {
			std::string line;
			std::getline(channelsAndResectionFile, line);
			std::stringstream ss(line);

			channelAndResection readChAndResect;
			readChAndResect.clear();
			ss >> readChAndResect.channelName;
			ss >> readChAndResect.soz;
			ss >> readChAndResect.electrodeType;
			ss >> readChAndResect.artifact;
			ss >> readChAndResect.resection;
			ss >> readChAndResect.engelOuctome;

			if (!readChAndResect.channelName.empty()) {
				std::string channelName = readChAndResect.channelName;
				channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '-'), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '_'), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '#'), channelName.end());
				readChAndResect.channelName = channelName;

				//m_selectedBipolarChanels.push_back(readChAndResect.channelName);
				//m_channelsHFO_Statistics.push_back(readChAndResect);

				int nrSelectedChanns = m_useMontage ? m_selectedBipolarChanels.size() : m_selectedUnipolarChanels.size();
				for (int selChIdx = 0; selChIdx < nrSelectedChanns; selChIdx++) {
					std::string selChannName = m_useMontage ? m_selectedBipolarChanels[selChIdx] : m_selectedUnipolarChanels[selChIdx];
					selChannName.erase(std::remove_if(selChannName.begin(), selChannName.end(), ::isspace), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '-'), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '_'), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '#'), selChannName.end());
					if (selChannName.compare(channelName) == 0) {	// is channel from resection file found within the channels from the EEG file?
						m_channelsHFO_Statistics.push_back(readChAndResect);
					}
				}
			}

		}

		m_selectedBipolarChanels.clear();
		for (int analysisChIdx = 0; analysisChIdx < m_channelsHFO_Statistics.size(); analysisChIdx++) {	// the selected EEG channels will be the same as the channles found in resection file
			m_selectedBipolarChanels.push_back(m_channelsHFO_Statistics[analysisChIdx].channelName);
		}

		channelsAndResectionFile.close();
	}

}

void HFO_Detector::getChannelsHFO_Statistics(std::string channelName, std::vector<EOI_Event> &readEvents) {

	channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
	channelName.erase(std::remove(channelName.begin(), channelName.end(), '-'), channelName.end());
	channelName.erase(std::remove(channelName.begin(), channelName.end(), '_'), channelName.end());
	channelName.erase(std::remove(channelName.begin(), channelName.end(), '#'), channelName.end());

	long nrDetections = readEvents.size();
	int nrAnalyzedChannels = m_channelsHFO_Statistics.size();

	for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {

		if (m_channelsHFO_Statistics[ch].channelName.compare(channelName) == 0) {	// find right channel

			for (long detIdx = 0; detIdx < nrDetections; ++detIdx) {				// read detection
				std::string label = readEvents[detIdx].description;
				bool foundFastRipple = label.find("ast") != std::string::npos;
				bool foundRipple = ((label.find("ipple") != std::string::npos) || (label.find("HFO") != std::string::npos)) && !foundFastRipple;
				bool foundRippleAndFastRipple = false;
				bool foundSpike = label.find("ike") != std::string::npos;
				bool foundSleepSpindle = label.find("indle") != std::string::npos;

				if (label.find("ipple") != std::string::npos) {
					label.erase(label.find("ipple"), 5);
					foundRippleAndFastRipple = label.find("ipple") != std::string::npos;
					if (foundRippleAndFastRipple) {
						foundRipple = false;
						foundFastRipple = false;
					}
				}

				//IES
				if (foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikes++;
				}

				//Ripples, FRs, HFO
				if (foundRipple) {
					m_channelsHFO_Statistics[ch].totalRipples++;
					m_channelsHFO_Statistics[ch].totalHFO++;
				} else if (foundFastRipple) {
					m_channelsHFO_Statistics[ch].totalFRs++;
					m_channelsHFO_Statistics[ch].totalHFO++;
				}
				else if (foundRippleAndFastRipple) {
					m_channelsHFO_Statistics[ch].totalRipples++;
					m_channelsHFO_Statistics[ch].totalFRs++;
					m_channelsHFO_Statistics[ch].totalHFO+=2;
				}

				// IES-Ripples, IES-FRs, IES-HFO
				if (foundRipple && foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeAndRipple++;
					m_channelsHFO_Statistics[ch].totalSpikeAndHFO++;
				}
				if (foundFastRipple && foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeAndFR++;
					m_channelsHFO_Statistics[ch].totalSpikeAndHFO++;
				}
				if (foundRippleAndFastRipple && foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeAndRipple++;
					m_channelsHFO_Statistics[ch].totalSpikeAndFR++;
					m_channelsHFO_Statistics[ch].totalSpikeAndHFO+=2;
				}

				// IESorRipples, IESorFRs, IESorHFO
				if (foundRipple || foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeOrRipple++;
					m_channelsHFO_Statistics[ch].totalSpikeOrHFO++;
				}
				if (foundFastRipple || foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeOrFR++;
					m_channelsHFO_Statistics[ch].totalSpikeOrHFO++;
				}
				if (foundRippleAndFastRipple || foundSpike) {
					m_channelsHFO_Statistics[ch].totalSpikeOrRipple++;
					m_channelsHFO_Statistics[ch].totalSpikeOrFR++;
					m_channelsHFO_Statistics[ch].totalSpikeOrHFO+=2;
				}

				m_channelsHFO_Statistics[ch].totalSpikeAndRipplePlusFRs = m_channelsHFO_Statistics[ch].totalSpikeAndRipple + m_channelsHFO_Statistics[ch].totalFRs;
			}

			m_channelsHFO_Statistics[ch].nrMinutes++;

			m_channelsHFO_Statistics[ch].rippleRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalRipples / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].frRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalFRs / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].hfoRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalHFO / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spikeRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikes / (double)m_channelsHFO_Statistics[ch].nrMinutes;

			m_channelsHFO_Statistics[ch].spikeAndRippleRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeAndRipple / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spikeAndFrRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeAndFR / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spikeAndHfoRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeAndHFO / (double)m_channelsHFO_Statistics[ch].nrMinutes;

			m_channelsHFO_Statistics[ch].spikeOrRippleRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeOrRipple / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spikeOrFrRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeOrFR / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spikeOrHfoRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeOrHFO / (double)m_channelsHFO_Statistics[ch].nrMinutes;

			m_channelsHFO_Statistics[ch].spikeAndRipplePlusFRs_RatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpikeAndRipplePlusFRs / (double)m_channelsHFO_Statistics[ch].nrMinutes;

			/*m_channelsHFO_Statistics[ch].spindleRippleRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpindleRipple / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spindleFrRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpindleFR / (double)m_channelsHFO_Statistics[ch].nrMinutes;
			m_channelsHFO_Statistics[ch].spindleHfoRatePerMinute = (double)m_channelsHFO_Statistics[ch].totalSpindleHFO / (double)m_channelsHFO_Statistics[ch].nrMinutes;*/

		}
	}
}

void HFO_Detector::getROCs(bool localizeResectionOrSOZ) {

	int nrAnalyzedChannels = m_channelsHFO_Statistics.size();

	//A
	for (int paramIdx = 1; paramIdx <= 19; ++paramIdx) {
		
		if (paramIdx >= 5 && paramIdx <= 12) {
			continue;
		}

		double thMaxValue, thStepSize;
		getMaxThAndStepSize_ROC(paramIdx, thMaxValue, thStepSize);
		std::vector<double> thVec;
		for (double th = 0; th <= thMaxValue; th += thStepSize) {
			thVec.push_back(th);
		}
		long vecSize = thVec.size();
		std::vector<double> sensVec(vecSize, 0), fprVec(vecSize, 0), precVec(vecSize, 0), mccVector(vecSize, 0), kappaVector(vecSize, 0);

		std::vector<double> paramValsPerChann;
		std::vector<bool> resectionPerChannel;

		double auc = 0;
		std::string paramName;
		//Get Parameters Name
		if (paramIdx == 1) paramName = "rippleRatePerMinute";
		if (paramIdx == 2) paramName = "frRatePerMinute";
		if (paramIdx == 3) paramName = "hfoRatePerMinute";
		if (paramIdx == 4) paramName = "spikeRatePerMinute";

		if (paramIdx == 5) paramName = "maxRippleRatePerMin";
		if (paramIdx == 6) paramName = "maxFRRatePerMin";
		if (paramIdx == 7) paramName = "maxHFO_RatePerMinute";
		if (paramIdx == 8) paramName = "maxSpikeRatePerMinute";

		if (paramIdx == 9) paramName = "totalRipples";
		if (paramIdx == 10) paramName = "totalFRs";
		if (paramIdx == 11) paramName = "totalHFO";
		if (paramIdx == 12) paramName = "totalSpikes";

		if (paramIdx == 13) paramName = "spikeOrRippleRatePerMinute";
		if (paramIdx == 14) paramName = "spikeOrFrRatePerMinute";
		if (paramIdx == 15) paramName = "spikeOrHfoRatePerMinute";

		if (paramIdx == 16) paramName = "spikeAndRippleRatePerMinute";
		if (paramIdx == 17) paramName = "spikeAndFrRatePerMinute";
		if (paramIdx == 18) paramName = "spikeAndHfoRatePerMinute";
		if (paramIdx == 19) paramName = "spikeAndRipplePlusFRs_RatePerMinute";

		if (paramIdx == 20) paramName = "spindleRippleRatePerMinute";
		if (paramIdx == 21) paramName = "spindleFrRatePerMinute";
		if (paramIdx == 22) paramName = "spindleHfoRatePerMinute";

		//Get Parameters value per channel
		for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
			bool resectedChannel = localizeResectionOrSOZ ? m_channelsHFO_Statistics[ch].resection > 0 : m_channelsHFO_Statistics[ch].soz > 0;
			double parameterVal = 0;

			if (paramIdx == 1) parameterVal = m_channelsHFO_Statistics[ch].rippleRatePerMinute;
			if (paramIdx == 2) parameterVal = m_channelsHFO_Statistics[ch].frRatePerMinute;
			if (paramIdx == 3) parameterVal = m_channelsHFO_Statistics[ch].hfoRatePerMinute;
			if (paramIdx == 4) parameterVal = m_channelsHFO_Statistics[ch].spikeRatePerMinute;

			if (paramIdx == 5) parameterVal = m_channelsHFO_Statistics[ch].maxRippleRatePerMin;
			if (paramIdx == 6) parameterVal = m_channelsHFO_Statistics[ch].maxFRRatePerMin;
			if (paramIdx == 7) parameterVal = m_channelsHFO_Statistics[ch].maxHFO_RatePerMinute;
			if (paramIdx == 8) parameterVal = m_channelsHFO_Statistics[ch].maxSpikeRatePerMinute;

			if (paramIdx == 9) parameterVal = m_channelsHFO_Statistics[ch].totalRipples;
			if (paramIdx == 10) parameterVal = m_channelsHFO_Statistics[ch].totalFRs;
			if (paramIdx == 11) parameterVal = m_channelsHFO_Statistics[ch].totalHFO;
			if (paramIdx == 12) parameterVal = m_channelsHFO_Statistics[ch].totalSpikes;

			if (paramIdx == 13) parameterVal = m_channelsHFO_Statistics[ch].spikeOrRippleRatePerMinute;
			if (paramIdx == 14) parameterVal = m_channelsHFO_Statistics[ch].spikeOrFrRatePerMinute;
			if (paramIdx == 15) parameterVal = m_channelsHFO_Statistics[ch].spikeOrHfoRatePerMinute;

			if (paramIdx == 16) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRippleRatePerMinute;
			if (paramIdx == 17) parameterVal = m_channelsHFO_Statistics[ch].spikeAndFrRatePerMinute;
			if (paramIdx == 18) parameterVal = m_channelsHFO_Statistics[ch].spikeAndHfoRatePerMinute;
			if (paramIdx == 19) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRipplePlusFRs_RatePerMinute;

			if (paramIdx == 20) parameterVal = m_channelsHFO_Statistics[ch].spindleRippleRatePerMinute;
			if (paramIdx == 21) parameterVal = m_channelsHFO_Statistics[ch].spindleFrRatePerMinute;
			if (paramIdx == 22) parameterVal = m_channelsHFO_Statistics[ch].spindleHfoRatePerMinute;

			resectionPerChannel.push_back(resectedChannel);
			paramValsPerChann.push_back(parameterVal);
		}

#pragma omp parallel for
		for (long thIdx = 0; thIdx < thVec.size(); thIdx++) { 	// Max 125 HFO (FR) per second, i.e. 7500 HFO per Minute
			double th = thVec[thIdx];
			double stepSens = 0, stepSpec = 0, stepFPR, stepPrec = 0;
			double tp = 0, fp = 0, tn = 0, fn = 0;
			
			for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
				bool resectedChannel = resectionPerChannel[ch];
				double parameterVal = paramValsPerChann[ch];

				if (parameterVal >= th && resectedChannel) {
					tp++;
				}
				if (parameterVal >= th && !resectedChannel) {
					fp++;
				}
				if (parameterVal < th && !resectedChannel) {
					tn++;
				}
				if (parameterVal < th && resectedChannel) {
					fn++;
				}
			}
			stepSens = (tp + fn) > 0 ? tp / (tp + fn) : 0;
			stepSpec = (tn + fp) > 0 ? tn / (tn + fp) : 0;
			stepFPR = 1 - stepSpec;
			stepPrec = (tp + fp) > 0 ? tp / (tp + fp) : 0;

			// Matthews Correlation Coefficient
			double matthewsCC_num = (tp * tn - fp * fn);
			double matthewsCC_denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
			double stepMCC = matthewsCC_denom > 0 ? matthewsCC_num / matthewsCC_denom : 0;

			//Kappa
			double sum = tp + fn + fp + tn;
			double p_o = (tp + tn) / sum;
			double marginalA = ((tp + fn)*(tp + fp)) / sum;
			double marginalB = ((fp + tn)*(fn + tn)) / sum;
			double p_e = (marginalA + marginalB) / sum;
			double stepKappa = (p_o - p_e) / (1 - p_e);

			sensVec[thIdx] = stepSens * 100;
			fprVec[thIdx] = stepFPR * 100;
			precVec[thIdx] = stepPrec * 100;
			mccVector[thIdx] = stepMCC * 100;
			kappaVector[thIdx] = stepKappa * 100;
		}

		double maxMCC = 0, maxMCC_Sens = 0, maxMCC_Spec = 0, maxMCC_Prec = 0, maxMCC_Th = 0, maxKappa = 0;
		for (int i = 0; i < sensVec.size() - 1; ++i) {
			auc += abs(fprVec[i] - fprVec[i + 1]) * ((sensVec[i] + sensVec[i + 1]) / 2);

			if (mccVector[i] > maxMCC) {
				maxMCC = mccVector[i];
				maxMCC_Sens = sensVec[i];
				maxMCC_Spec = 1 - fprVec[i];
				maxMCC_Prec = precVec[i];
				maxMCC_Th = thVec[i];
				maxKappa = kappaVector[i];
			}
		}
		auc /= 100;
		//writeROC_File(localizeResectionOrSOZ, paramName, m_channelsHFO_Statistics[0].nrMinutes, thVec, fprVec, sensVec, precVec, mccVector, kappaVector, auc);
		//writeAUC_File(localizeResectionOrSOZ, m_channelsHFO_Statistics[0].nrMinutes, paramName, auc, maxMCC, maxMCC_Sens, maxMCC_Spec, maxMCC_Prec, maxMCC_Th, maxKappa);
		writeCommonAUC_File(localizeResectionOrSOZ, m_channelsHFO_Statistics[0].nrMinutes, paramName, auc, maxMCC, maxMCC_Sens, maxMCC_Spec, maxMCC_Prec, maxMCC_Th, maxKappa);
		writeInterPatNormParamVals_File(localizeResectionOrSOZ, m_channelsHFO_Statistics[0].nrMinutes, paramName, resectionPerChannel, paramValsPerChann);

	}
}

/*void HFO_Detector::getPercentageROCs(bool localizeResectionOrSOZ) {

	int nrAnalyzedChannels = m_channelsHFO_Statistics.size();

	//A
	for (int paramIdx = 1; paramIdx <= 19; ++paramIdx) {
		std::vector<double> thVec, sensVec, fprVec, precVec, mccVector, kappaVector;
		double auc = 0;
		std::string paramName;

		if (paramIdx == 1) paramName = "rippleRatePerMinute";
		if (paramIdx == 2) paramName = "frRatePerMinute";
		if (paramIdx == 3) paramName = "hfoRatePerMinute";
		if (paramIdx == 4) paramName = "spikeRatePerMinute";

		if (paramIdx == 5) paramName = "maxRippleRatePerMin";
		if (paramIdx == 6) paramName = "maxFRRatePerMin";
		if (paramIdx == 7) paramName = "maxHFO_RatePerMinute";
		if (paramIdx == 8) paramName = "maxSpikeRatePerMinute";

		if (paramIdx == 9) paramName = "totalRipples";
		if (paramIdx == 10) paramName = "totalFRs";
		if (paramIdx == 11) paramName = "totalHFO";
		if (paramIdx == 12) paramName = "totalSpikes";

		if (paramIdx == 13) paramName = "spikeOrRippleRatePerMinute";
		if (paramIdx == 14) paramName = "spikeOrFrRatePerMinute";
		if (paramIdx == 15) paramName = "spikeOrHfoRatePerMinute";

		if (paramIdx == 16) paramName = "spikeAndRippleRatePerMinute";
		if (paramIdx == 17) paramName = "spikeAndFrRatePerMinute";
		if (paramIdx == 18) paramName = "spikeAndHfoRatePerMinute";
		if (paramIdx == 19) paramName = "spikeAndRipplePlusFRs_RatePerMinute";

		if (paramIdx == 20) paramName = "spindleRippleRatePerMinute";
		if (paramIdx == 21) paramName = "spindleFrRatePerMinute";
		if (paramIdx == 22) paramName = "spindleHfoRatePerMinute";

		paramName += "_Percentage_";

		double maxAcrossChanns = 0;
		for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
			double parameterVal = 0;
			if (paramIdx == 1) parameterVal = m_channelsHFO_Statistics[ch].rippleRatePerMinute;
			if (paramIdx == 2) parameterVal = m_channelsHFO_Statistics[ch].frRatePerMinute;
			if (paramIdx == 3) parameterVal = m_channelsHFO_Statistics[ch].hfoRatePerMinute;
			if (paramIdx == 4) parameterVal = m_channelsHFO_Statistics[ch].spikeRatePerMinute;

			if (paramIdx == 5) parameterVal = m_channelsHFO_Statistics[ch].maxRippleRatePerMin;
			if (paramIdx == 6) parameterVal = m_channelsHFO_Statistics[ch].maxFRRatePerMin;
			if (paramIdx == 7) parameterVal = m_channelsHFO_Statistics[ch].maxHFO_RatePerMinute;
			if (paramIdx == 8) parameterVal = m_channelsHFO_Statistics[ch].maxSpikeRatePerMinute;

			if (paramIdx == 9) parameterVal = m_channelsHFO_Statistics[ch].totalRipples;
			if (paramIdx == 10) parameterVal = m_channelsHFO_Statistics[ch].totalFRs;
			if (paramIdx == 11) parameterVal = m_channelsHFO_Statistics[ch].totalHFO;
			if (paramIdx == 12) parameterVal = m_channelsHFO_Statistics[ch].totalSpikes;

			if (paramIdx == 13) parameterVal = m_channelsHFO_Statistics[ch].spikeOrRippleRatePerMinute;
			if (paramIdx == 14) parameterVal = m_channelsHFO_Statistics[ch].spikeOrFrRatePerMinute;
			if (paramIdx == 15) parameterVal = m_channelsHFO_Statistics[ch].spikeOrHfoRatePerMinute;

			if (paramIdx == 16) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRippleRatePerMinute;
			if (paramIdx == 17) parameterVal = m_channelsHFO_Statistics[ch].spikeAndFrRatePerMinute;
			if (paramIdx == 18) parameterVal = m_channelsHFO_Statistics[ch].spikeAndHfoRatePerMinute;
			if (paramIdx == 19) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRipplePlusFRs_RatePerMinute;


			if (paramIdx == 20) parameterVal = m_channelsHFO_Statistics[ch].spindleRippleRatePerMinute;
			if (paramIdx == 21) parameterVal = m_channelsHFO_Statistics[ch].spindleFrRatePerMinute;
			if (paramIdx == 22) parameterVal = m_channelsHFO_Statistics[ch].spindleHfoRatePerMinute;

			if (parameterVal > maxAcrossChanns) {
				maxAcrossChanns = parameterVal;
			}
		}

		for (double th = 0; th <= 1; th += 0.01) { 	// Max 125 HFO (FR) per second, i.e. 7500 HFO per Minute
			double stepSens = 0, stepSpec = 0, stepFPR, stepPrec = 0;
			double tp = 0, fp = 0, tn = 0, fn = 0;

			for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
				bool resectedChannel = localizeResectionOrSOZ ? m_channelsHFO_Statistics[ch].resection > 0 : m_channelsHFO_Statistics[ch].soz > 0;

				double parameterVal = 0;
				if (paramIdx == 1) parameterVal = m_channelsHFO_Statistics[ch].rippleRatePerMinute / maxAcrossChanns;
				if (paramIdx == 2) parameterVal = m_channelsHFO_Statistics[ch].frRatePerMinute / maxAcrossChanns;
				if (paramIdx == 3) parameterVal = m_channelsHFO_Statistics[ch].hfoRatePerMinute / maxAcrossChanns;
				if (paramIdx == 4) parameterVal = m_channelsHFO_Statistics[ch].spikeRatePerMinute / maxAcrossChanns;

				if (paramIdx == 5) parameterVal = m_channelsHFO_Statistics[ch].maxRippleRatePerMin / maxAcrossChanns;
				if (paramIdx == 6) parameterVal = m_channelsHFO_Statistics[ch].maxFRRatePerMin / maxAcrossChanns;
				if (paramIdx == 7) parameterVal = m_channelsHFO_Statistics[ch].maxHFO_RatePerMinute / maxAcrossChanns;
				if (paramIdx == 8) parameterVal = m_channelsHFO_Statistics[ch].maxSpikeRatePerMinute / maxAcrossChanns;

				if (paramIdx == 9) parameterVal = m_channelsHFO_Statistics[ch].totalRipples / maxAcrossChanns;
				if (paramIdx == 10) parameterVal = m_channelsHFO_Statistics[ch].totalFRs / maxAcrossChanns;
				if (paramIdx == 11) parameterVal = m_channelsHFO_Statistics[ch].totalHFO / maxAcrossChanns;
				if (paramIdx == 12) parameterVal = m_channelsHFO_Statistics[ch].totalSpikes / maxAcrossChanns;

				if (paramIdx == 13) parameterVal = m_channelsHFO_Statistics[ch].spikeOrRippleRatePerMinute / maxAcrossChanns;
				if (paramIdx == 14) parameterVal = m_channelsHFO_Statistics[ch].spikeOrFrRatePerMinute / maxAcrossChanns;
				if (paramIdx == 15) parameterVal = m_channelsHFO_Statistics[ch].spikeOrHfoRatePerMinute / maxAcrossChanns;

				if (paramIdx == 16) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRippleRatePerMinute / maxAcrossChanns;
				if (paramIdx == 17) parameterVal = m_channelsHFO_Statistics[ch].spikeAndFrRatePerMinute / maxAcrossChanns;
				if (paramIdx == 18) parameterVal = m_channelsHFO_Statistics[ch].spikeAndHfoRatePerMinute / maxAcrossChanns;
				if (paramIdx == 19) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRipplePlusFRs_RatePerMinute / maxAcrossChanns;

				if (paramIdx == 20) parameterVal = m_channelsHFO_Statistics[ch].spindleRippleRatePerMinute / maxAcrossChanns;
				if (paramIdx == 21) parameterVal = m_channelsHFO_Statistics[ch].spindleFrRatePerMinute / maxAcrossChanns;
				if (paramIdx == 22) parameterVal = m_channelsHFO_Statistics[ch].spindleHfoRatePerMinute / maxAcrossChanns;


				if (parameterVal >= th && resectedChannel) {
					tp++;
				}
				if (parameterVal >= th && !resectedChannel) {
					fp++;
				}
				if (parameterVal < th && !resectedChannel) {
					tn++;
				}
				if (parameterVal < th && resectedChannel) {
					fn++;
				}
			}
			stepSens = (tp + fn) > 0 ? tp / (tp + fn) : 0;
			stepSpec = (tn + fp) > 0 ? tn / (tn + fp) : 0;
			stepFPR = 1 - stepSpec;
			stepPrec = (tp + fp) > 0 ? tp / (tp + fp) : 0;

			// Matthews Correlation Coefficient
			double matthewsCC_num = (tp * tn - fp * fn);
			double matthewsCC_denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
			double stepMCC = matthewsCC_denom > 0 ? matthewsCC_num / matthewsCC_denom : 0;

			//Kappa
			double sum = tp + fn + fp + tn;
			double p_o = (tp + tn) / sum;
			double marginalA = ((tp + fn)*(tp + fp)) / sum;
			double marginalB = ((fp + tn)*(fn + tn)) / sum;
			double p_e = (marginalA + marginalB) / sum;
			double stepKappa = (p_o - p_e) / (1 - p_e);

			thVec.push_back(th);
			sensVec.push_back(stepSens*100);
			fprVec.push_back(stepFPR * 100);
			precVec.push_back(stepPrec * 100);
			mccVector.push_back(stepMCC * 100);
			kappaVector.push_back(stepKappa * 100);
		}
		double maxMCC = 0, maxMCC_Sens = 0, maxMCC_Spec = 0, maxMCC_Prec = 0, maxMCC_Th = 0, maxKappa = 0;
		for (int i = 0; i < sensVec.size() - 1; ++i) {
			auc += abs(fprVec[i] - fprVec[i + 1]) * ((sensVec[i] + sensVec[i + 1]) / 2);

			if (mccVector[i] > maxMCC) {
				maxMCC = mccVector[i];
				maxMCC_Sens = sensVec[i];
				maxMCC_Spec = 1 - fprVec[i];
				maxMCC_Prec = precVec[i];
				maxMCC_Th = thVec[i];
				maxKappa = kappaVector[i];
			}
		}

		writeROC_File(localizeResectionOrSOZ, paramName, m_channelsHFO_Statistics[0].nrMinutes, thVec, fprVec, sensVec, precVec, mccVector, kappaVector, auc);
		writeAUC_File(localizeResectionOrSOZ, m_channelsHFO_Statistics[0].nrMinutes, paramName, auc, maxMCC, maxMCC_Sens, maxMCC_Spec, maxMCC_Prec, maxMCC_Th, maxKappa);
	}
}*/

void HFO_Detector::writeROC_File(bool localizeResectionOrSOZ, std::string parameterName, int minute, std::vector<double> thVec, std::vector<double> fprVec, std::vector<double> sensVec, std::vector<double> precVec, std::vector<double> mccVector, std::vector<double> kappaVector, double areUnderCurve) {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName + "/";
	_mkdir(patPath.c_str());
	std::string rocPath = patPath + "/" + "ROC/";
	_mkdir(rocPath.c_str());
	rocPath += localizeResectionOrSOZ ? "Resection/" : "SOZ/";
	_mkdir(rocPath.c_str());

	std::string filename = rocPath + "ROC_" + std::to_string(minute) + "_" + parameterName + "_AUC=" + std::to_string(areUnderCurve) +  ".txt";

	std::ofstream vectorsFile;
	std::ifstream infile(filename);
	if (!infile.good()) {
		//Header
		vectorsFile.open(filename);
		vectorsFile << "Threshold" << "\t" << "FPR" << "\t" << "Sensitivity" << "\t" << "Precision" << "\t" << "MCC" << "\t" << "Kappa" << "\t" << areUnderCurve << "\n";
	}
	else {
		vectorsFile.open(filename, std::ios::app);    // open file for appending
	}

	for (long i = 0; i < fprVec.size(); ++i) {
		vectorsFile << thVec[i] << "\t" << fprVec[i] << "\t" << sensVec[i] << "\t" << precVec[i] << "\t" << mccVector[i] << "\t" << kappaVector[i] << "\n";
	}

	vectorsFile.close();
}

void HFO_Detector::writeAUC_File(bool localizeResectionOrSOZ, int minute, std::string parameterName, double areUnderCurve, double bestMCC, double maxMCC_Sens, double maxMCC_Spec, double maxMCC_Prec, double maxMCC_Th, double maxKappa) {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	std::string patPath = outputPath + "/" + m_patientName + "/";
	_mkdir(patPath.c_str());
	std::string rocPath = patPath + "/" + "ROC/";
	_mkdir(rocPath.c_str());
	rocPath += localizeResectionOrSOZ ? "Resection/" : "SOZ/";
	_mkdir(rocPath.c_str());

	std::string filename = rocPath + "AreasUnderTheCurve.txt";

	std::ofstream vectorsFile;
	std::ifstream infile(filename);
	if (!infile.good()) {
		//Header
		vectorsFile.open(filename);
		vectorsFile << "Minute" << "\t" << "parameterName" << "\t" << "areUnderCurve" << "\t" << "bestMCC" << "\t" << "maxMCC_Sens" << "\t" << "maxMCC_Spec" << "\t" << "maxMCC_Prec" << "\t" << "maxMCC_Th" << "\t" << "maxKappa" << "\n";
	}
	else {
		vectorsFile.open(filename, std::ios::app);    // open file for appending
	}

	vectorsFile << minute << "\t" << parameterName << "\t" << areUnderCurve << "\t" << bestMCC << "\t" << maxMCC_Sens << "\t" << maxMCC_Spec << "\t" << maxMCC_Prec << "\t" << maxMCC_Th << "\t" << maxKappa << "\n";;

	vectorsFile.close();
}

void HFO_Detector::writeInterPatNormParamVals_File(bool localizeResectionOrSOZ, int minute, std::string paramName, std::vector<bool> resectionPerChannel, std::vector<double> paramValsPerChann) {
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";

	std::string goldStandard = localizeResectionOrSOZ ? "Resection" : "SOZ";
	std::string filename = outputPath + "interPatNormParamVals_" + paramName + "_" + goldStandard + ".txt";

	std::ofstream vectorsFile;
	std::ifstream infile(filename);
	if (!infile.good()) {
		//Header
		vectorsFile.open(filename);
		vectorsFile << "patientName" << "\t" << "minute" << "\t" << paramName << "\t" << goldStandard << "\n";
	}
	else {
		vectorsFile.open(filename, std::ios::app);    // open file for appending
	}

	//Normalize
	vectorStats paramValsStats;
	paramValsStats = SignalProcessing::vectorStatistics(paramValsPerChann);
	for (int i = 0; i < paramValsPerChann.size(); ++i) {
		double normVal = (paramValsPerChann[i] - paramValsStats.mean) / paramValsStats.stdrdDev;
		vectorsFile << m_patientName << "\t" << minute << "\t" << normVal << "\t" << resectionPerChannel[i] << "\n";
	}
	vectorsFile.close();
}

void HFO_Detector::writeCommonAUC_File(bool localizeResectionOrSOZ, int minute, std::string parameterName, double areUnderCurve, double bestMCC, double maxMCC_Sens, double maxMCC_Spec, double maxMCC_Prec, double maxMCC_Th, double maxKappa) {
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	//std::string patPath = outputPath + "/" + m_patientName + "/";
	std::string goldStandard = localizeResectionOrSOZ ? "Resection" : "SOZ";
	std::string filename = outputPath + "interPatientAUC_" + parameterName + "_" + goldStandard + ".txt";

	std::ofstream commonAUC_File;
	std::ifstream infile(filename);
	if (!infile.good()) {
		//Header
		commonAUC_File.open(filename);
		commonAUC_File << "Parameter" << "\t" << "Minute" << "\t" << "PatientName" << "\t" << "areUnderCurve" << "\t" << "bestMCC" << "\t" << "maxMCC_Sens" << "\t" << "maxMCC_Spec" << "\t" << "maxMCC_Prec" << "\t" << "maxMCC_Th" << "\t" << "maxKappa" << "\n";
	}
	else {
		commonAUC_File.open(filename, std::ios::app);    // open file for appending
	}

	commonAUC_File << parameterName << "\t" << minute << "\t" << m_patientName << "\t" << areUnderCurve << "\t" << bestMCC << "\t" << maxMCC_Sens << "\t" << maxMCC_Spec << "\t" << maxMCC_Prec << "\t" << maxMCC_Th << "\t" << maxKappa << "\n";;

	commonAUC_File.close();
}

void HFO_Detector::writeDetectionsForSections() {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	std::string filename = outputPath + "detectionsForSectionsFile.txt";
	remove(filename.c_str());

	std::ofstream sectionsFile;
	sectionsFile.open(filename);
	sectionsFile << "nrMinutes" << "\t" << "ChannelName" << "\t" << "totalRipples" << "\t" << "totalFRs" << "\t" << "totalHFO" << "\t" << "totalSpikes" << "\t"
		<< "totalSpikeOrRipple" << "\t" << "totalSpikeOrFR" << "\t" << "totalSpikeOrHFO" << "\t"
		<< "totalSpikeAndRipple" << "\t" << "totalSpikeAndFR" << "\t" << "totalSpikeAndHFO" << "\t" << "totalSpikeAndRipplePlusFRs" << "\t"
		<< "totalSpindleRipple" << "\t" << "totalSpindleFR" << "\t" << "totalSpindleHFO" << "\n";

	int nrAnalyzedChannels = m_channelsHFO_Statistics.size();
	for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
		sectionsFile << m_channelsHFO_Statistics[ch].nrMinutes << "\t" << m_channelsHFO_Statistics[ch].channelName << "\t"
			<< m_channelsHFO_Statistics[ch].totalRipples << "\t" << m_channelsHFO_Statistics[ch].totalFRs << "\t" << m_channelsHFO_Statistics[ch].totalHFO << "\t" << m_channelsHFO_Statistics[ch].totalSpikes << "\t"
			<< m_channelsHFO_Statistics[ch].totalSpikeOrRipple << "\t" << m_channelsHFO_Statistics[ch].totalSpikeOrFR << "\t" << m_channelsHFO_Statistics[ch].totalSpikeOrHFO << "\t"
			<< m_channelsHFO_Statistics[ch].totalSpikeAndRipple << "\t" << m_channelsHFO_Statistics[ch].totalSpikeAndFR << "\t" << m_channelsHFO_Statistics[ch].totalSpikeAndHFO << "\t" << m_channelsHFO_Statistics[ch].totalSpikeAndRipplePlusFRs << "\t"
			<< m_channelsHFO_Statistics[ch].totalSpindleRipple << "\t" << m_channelsHFO_Statistics[ch].totalSpindleFR << "\t" << m_channelsHFO_Statistics[ch].totalSpindleHFO << "\n";
	}
	sectionsFile.close();
}

void HFO_Detector::readDetectionsForSections() {
	
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	std::string fn = outputPath + "detectionsForSectionsFile.txt";

	if (fileExists(fn.c_str())) {

		std::ifstream sectionsFile;
		sectionsFile.open(fn.c_str());
		std::string firstLine;
		std::getline(sectionsFile, firstLine);

		while (!sectionsFile.eof()) {
			std::string line;
			std::getline(sectionsFile, line);
			std::stringstream ss(line);
			channelAndResection readSavedChannelAndResectionData;

			ss >> readSavedChannelAndResectionData.nrMinutes;
			ss >> readSavedChannelAndResectionData.channelName;

			ss >> readSavedChannelAndResectionData.totalRipples;
			ss >> readSavedChannelAndResectionData.totalFRs;
			ss >> readSavedChannelAndResectionData.totalHFO;
			ss >> readSavedChannelAndResectionData.totalSpikes;

			ss >> readSavedChannelAndResectionData.totalSpikeOrRipple;
			ss >> readSavedChannelAndResectionData.totalSpikeOrFR;
			ss >> readSavedChannelAndResectionData.totalSpikeOrHFO;

			ss >> readSavedChannelAndResectionData.totalSpikeAndRipple;
			ss >> readSavedChannelAndResectionData.totalSpikeAndFR;
			ss >> readSavedChannelAndResectionData.totalSpikeAndHFO;
			ss >> readSavedChannelAndResectionData.totalSpikeAndRipplePlusFRs;

			ss >> readSavedChannelAndResectionData.totalSpindleRipple;
			ss >> readSavedChannelAndResectionData.totalSpindleFR;
			ss >> readSavedChannelAndResectionData.totalSpindleHFO;

			if (!readSavedChannelAndResectionData.channelName.empty()) {
				std::string channelName = readSavedChannelAndResectionData.channelName;
				channelName.erase(std::remove_if(channelName.begin(), channelName.end(), ::isspace), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '-'), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '_'), channelName.end());
				channelName.erase(std::remove(channelName.begin(), channelName.end(), '#'), channelName.end());

				int nrAnalyzedChannels = m_channelsHFO_Statistics.size();
				for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
					std::string selChannName = m_channelsHFO_Statistics[ch].channelName;
					selChannName.erase(std::remove_if(selChannName.begin(), selChannName.end(), ::isspace), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '-'), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '_'), selChannName.end());
					selChannName.erase(std::remove(selChannName.begin(), selChannName.end(), '#'), selChannName.end());
					if (selChannName.compare(channelName) == 0) {	// is channel from resection file found within the channels from the EEG file?
						m_channelsHFO_Statistics[ch].nrMinutes = readSavedChannelAndResectionData.nrMinutes;

						m_channelsHFO_Statistics[ch].totalRipples = readSavedChannelAndResectionData.totalRipples;
						m_channelsHFO_Statistics[ch].totalFRs = readSavedChannelAndResectionData.totalFRs;
						m_channelsHFO_Statistics[ch].totalHFO = readSavedChannelAndResectionData.totalHFO;
						m_channelsHFO_Statistics[ch].totalSpikes = readSavedChannelAndResectionData.totalSpikes;

						m_channelsHFO_Statistics[ch].totalSpikeOrRipple = readSavedChannelAndResectionData.totalSpikeOrRipple;
						m_channelsHFO_Statistics[ch].totalSpikeOrFR = readSavedChannelAndResectionData.totalSpikeOrFR;
						m_channelsHFO_Statistics[ch].totalSpikeOrHFO = readSavedChannelAndResectionData.totalSpikeOrHFO;

						m_channelsHFO_Statistics[ch].totalSpikeAndRipple = readSavedChannelAndResectionData.totalSpikeAndRipple;
						m_channelsHFO_Statistics[ch].totalSpikeAndFR = readSavedChannelAndResectionData.totalSpikeAndFR;
						m_channelsHFO_Statistics[ch].totalSpikeAndHFO = readSavedChannelAndResectionData.totalSpikeAndHFO;
						m_channelsHFO_Statistics[ch].totalSpikeAndRipplePlusFRs = readSavedChannelAndResectionData.totalSpikeAndRipplePlusFRs;

						m_channelsHFO_Statistics[ch].totalSpindleRipple = readSavedChannelAndResectionData.totalSpindleRipple;
						m_channelsHFO_Statistics[ch].totalSpindleFR = readSavedChannelAndResectionData.totalSpindleFR;
						m_channelsHFO_Statistics[ch].totalSpindleHFO = readSavedChannelAndResectionData.totalSpindleHFO;
					}
				}
			}
		}
		sectionsFile.close();
	}
}


void HFO_Detector::getMaxThAndStepSize_ROC(int paramIdx, double &maxTh, double &stepSize) {
	int nrAnalyzedChannels = m_channelsHFO_Statistics.size();
	
	maxTh = 0;
	std::vector<double> paramVals;
	for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
		double parameterVal = 0;
		if (paramIdx == 1) parameterVal = m_channelsHFO_Statistics[ch].rippleRatePerMinute;
		if (paramIdx == 2) parameterVal = m_channelsHFO_Statistics[ch].frRatePerMinute;
		if (paramIdx == 3) parameterVal = m_channelsHFO_Statistics[ch].hfoRatePerMinute;
		if (paramIdx == 4) parameterVal = m_channelsHFO_Statistics[ch].spikeRatePerMinute;

		if (paramIdx == 5) parameterVal = m_channelsHFO_Statistics[ch].maxRippleRatePerMin;
		if (paramIdx == 6) parameterVal = m_channelsHFO_Statistics[ch].maxFRRatePerMin;
		if (paramIdx == 7) parameterVal = m_channelsHFO_Statistics[ch].maxHFO_RatePerMinute;
		if (paramIdx == 8) parameterVal = m_channelsHFO_Statistics[ch].maxSpikeRatePerMinute;

		if (paramIdx == 9) parameterVal = m_channelsHFO_Statistics[ch].totalRipples;
		if (paramIdx == 10) parameterVal = m_channelsHFO_Statistics[ch].totalFRs;
		if (paramIdx == 11) parameterVal = m_channelsHFO_Statistics[ch].totalHFO;
		if (paramIdx == 12) parameterVal = m_channelsHFO_Statistics[ch].totalSpikes;

		if (paramIdx == 13) parameterVal = m_channelsHFO_Statistics[ch].spikeOrRippleRatePerMinute;
		if (paramIdx == 14) parameterVal = m_channelsHFO_Statistics[ch].spikeOrFrRatePerMinute;
		if (paramIdx == 15) parameterVal = m_channelsHFO_Statistics[ch].spikeOrHfoRatePerMinute;

		if (paramIdx == 16) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRippleRatePerMinute;
		if (paramIdx == 17) parameterVal = m_channelsHFO_Statistics[ch].spikeAndFrRatePerMinute;
		if (paramIdx == 18) parameterVal = m_channelsHFO_Statistics[ch].spikeAndHfoRatePerMinute;
		if (paramIdx == 19) parameterVal = m_channelsHFO_Statistics[ch].spikeAndRipplePlusFRs_RatePerMinute;

		if (paramIdx == 20) parameterVal = m_channelsHFO_Statistics[ch].spindleRippleRatePerMinute;
		if (paramIdx == 21) parameterVal = m_channelsHFO_Statistics[ch].spindleFrRatePerMinute;
		if (paramIdx == 22) parameterVal = m_channelsHFO_Statistics[ch].spindleHfoRatePerMinute;

		paramVals.push_back(parameterVal);

		if (parameterVal > maxTh) {
			maxTh = parameterVal;
		}
	}

	std::sort(paramVals.begin(), paramVals.end());
	double minDiff = 1e10;
	for (int i = 0; i < paramVals.size()-1; ++i) {
		double diff = abs(paramVals[i] - paramVals[i + 1]);
		if (diff > 0.0 && diff < minDiff) {
			minDiff = diff;
		}
	}
	stepSize = minDiff;

	maxTh += 2* stepSize;
	stepSize /= 2;
}

void HFO_Detector::generatePerChannelAnnotations() {

	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	std::string patPath = outputPath + m_patientName;
	std::string detectionsPath = patPath + "/DetectionFiles/";
	std::vector<std::string> detTypes = { "RippleDetections", "FastRippleDetections" , "HFODetections" , "SpikeDetections" };


	for (int detTypeIdx = 0; detTypeIdx < detTypes.size(); detTypeIdx++) {

		std::unordered_map<std::string, std::vector<EOI_Event>> detectionsPerChannel;
		std::string readFilename = detectionsPath + m_patientName + "_" + detTypes[detTypeIdx] + ".txt";

		std::ifstream eoiDetections;
		eoiDetections.open(readFilename.c_str());
		eoiDetections.precision(32);
		while (!eoiDetections.eof()) {
			std::string line;
			std::getline(eoiDetections, line);
			std::stringstream ss(line);

			EOI_Event eoiDetection;
			std::string data;

			ss >> data; //:d
			if (data.compare("d:") == 0) {
				ss >> eoiDetection.description;
				ss >> eoiDetection.channelName;
				ss >> eoiDetection.channelNr;
				ss >> eoiDetection.startTime;
				ss >> eoiDetection.endTime;
				eoiDetection.duration = eoiDetection.endTime - eoiDetection.startTime;

				std::size_t found = eoiDetection.channelName.find("EEG ");
				if (found != std::string::npos) {
					eoiDetection.channelName.erase(found, found + 4);
				}

				found = eoiDetection.channelName.find("EEG");
				if (found != std::string::npos) {
					eoiDetection.channelName.erase(found, found + 3);
				}

				if (!eoiDetection.channelName.empty()) {
					if (detectionsPerChannel.find(eoiDetection.channelName) == detectionsPerChannel.end()) {
						detectionsPerChannel.emplace(eoiDetection.channelName, *new std::vector<EOI_Event>);
					}
					detectionsPerChannel.at(eoiDetection.channelName).push_back(eoiDetection);
				}
			}
		}
		eoiDetections.close();

		//write detections file for each channel
		for (auto& x : detectionsPerChannel) {
			std::string channName = x.first;
			writePerChannelDetectionsFile(detectionsPath, detTypes[detTypeIdx], x.first, x.second);

		}
	}
}

bool HFO_Detector::writePerChannelDetectionsFile(std::string directory, std::string eventType, std::string channName, std::vector<EOI_Event> &detections) {

	std::string detectionTypeDirectory = directory + eventType;
	_mkdir(detectionTypeDirectory.c_str());

	std::string filename = detectionTypeDirectory + "/" + m_patientName + "_" + channName + "_" + eventType + ".txt";
	
	double nrDetections = detections.size();
	std::vector<bool> touchedDetects(nrDetections, false);


	std::ofstream detectsFile;
	detectsFile.open(filename);
	//Header
	detectsFile << "SignalDuration:\t" << 0 << "\n";
	detectsFile << "Data\t" << "Description\t" << "ChannelName\t" << "Channel\t" << "StartTime(s)\t" << "EndTime(s)\t" << "StartTime(day)\t" << "RealEndTime(day)\t" << "Duration(s)\t" << "Match\n";

	for (long di = 0; di < nrDetections; di++) {

		detectsFile << "d:\t" <<
			detections[di].description << "\t" <<
			detections[di].channelName << "\t" <<
			detections[di].channelNr << "\t" <<
			detections[di].startTime << "\t" <<
			detections[di].endTime << "\t" <<
			0 << ":" << 0 << ":" << 0 << "\t" <<
			0 << ":" << 0 << ":" << 0 << "\t" <<
			detections[di].duration << "\t" <<
			0 << "\n";

	}

	detectsFile.close();
	return true;
}

void HFO_Detector::setBiomarkerOccChannels(void) {
	int nrAnalyzedChannels = m_useMontage ? m_montageLabels.size() : m_unipolarLabels.size();
	m_biomarkerOccPerChann.clear();
	m_biomarkerOccPerChann.resize(nrAnalyzedChannels);
	for (int ch = 0; ch < nrAnalyzedChannels; ++ch) {
		m_biomarkerOccPerChann[ch].clear();
		std::string channName = m_useMontage ? m_montageLabels[ch].montageName : m_unipolarLabels[ch].contactName;
		m_biomarkerOccPerChann[ch].channelName = channName;
	}
}

void HFO_Detector::getBiomarkerOcurrenceRates(int channIdx, std::vector<EOI_Event> &readEvents) {

	long nrDetections = readEvents.size();

	for (long detIdx = 0; detIdx < nrDetections; ++detIdx) {				// read detection, all of these detection correspond to the channel signaled by channIdx
		std::string label = readEvents[detIdx].description;
		bool foundFastRipple = label.find("ast") != std::string::npos;
		bool foundRipple = ((label.find("ipple") != std::string::npos) || (label.find("HFO") != std::string::npos)) && !foundFastRipple;
		bool foundRippleAndFastRipple = false;
		bool foundSpike = label.find("ike") != std::string::npos;
		bool foundSleepSpindle = label.find("indle") != std::string::npos;

		if (label.find("ipple") != std::string::npos) {
			label.erase(label.find("ipple"), 5);
			foundRippleAndFastRipple = label.find("ipple") != std::string::npos;
			if (foundRippleAndFastRipple) {
				foundRipple = false;
				foundFastRipple = false;
			}
		}

		//IES
		if (foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikes++;
		}

		//Ripples, FRs, HFO
		if (foundRipple) {
			m_biomarkerOccPerChann[channIdx].totalRipples++;
			m_biomarkerOccPerChann[channIdx].totalHFO++;
		}
		else if (foundFastRipple) {
			m_biomarkerOccPerChann[channIdx].totalFRs++;
			m_biomarkerOccPerChann[channIdx].totalHFO++;
		}
		else if (foundRippleAndFastRipple) {
			m_biomarkerOccPerChann[channIdx].totalRipples++;
			m_biomarkerOccPerChann[channIdx].totalFRs++;
			m_biomarkerOccPerChann[channIdx].totalHFO += 2;
		}

		// IES-Ripples, IES-FRs, IES-HFO
		if (foundRipple && foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeAndRipple++;
			m_biomarkerOccPerChann[channIdx].totalSpikeAndHFO++;
		}
		if (foundFastRipple && foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeAndFR++;
			m_biomarkerOccPerChann[channIdx].totalSpikeAndHFO++;
		}
		if (foundRippleAndFastRipple && foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeAndRipple++;
			m_biomarkerOccPerChann[channIdx].totalSpikeAndFR++;
			m_biomarkerOccPerChann[channIdx].totalSpikeAndHFO += 2;
		}

		// IESorRipples, IESorFRs, IESorHFO
		if (foundRipple || foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeOrRipple++;
			m_biomarkerOccPerChann[channIdx].totalSpikeOrHFO++;
		}
		if (foundFastRipple || foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeOrFR++;
			m_biomarkerOccPerChann[channIdx].totalSpikeOrHFO++;
		}
		if (foundRippleAndFastRipple || foundSpike) {
			m_biomarkerOccPerChann[channIdx].totalSpikeOrRipple++;
			m_biomarkerOccPerChann[channIdx].totalSpikeOrFR++;
			m_biomarkerOccPerChann[channIdx].totalSpikeOrHFO += 2;
		}

		m_biomarkerOccPerChann[channIdx].totalSpikeAndRipplePlusFRs = m_biomarkerOccPerChann[channIdx].totalSpikeAndRipple + m_biomarkerOccPerChann[channIdx].totalFRs;
	}

	m_biomarkerOccPerChann[channIdx].nrMinutes++;

	m_biomarkerOccPerChann[channIdx].rippleRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalRipples / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].frRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalFRs / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].hfoRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalHFO / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spikeRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikes / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;

	m_biomarkerOccPerChann[channIdx].spikeAndRippleRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeAndRipple / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spikeAndFrRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeAndFR / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spikeAndHfoRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeAndHFO / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;

	m_biomarkerOccPerChann[channIdx].spikeOrRippleRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeOrRipple / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spikeOrFrRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeOrFR / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spikeOrHfoRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeOrHFO / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;

	m_biomarkerOccPerChann[channIdx].spikeAndRipplePlusFRs_RatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpikeAndRipplePlusFRs / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;

	/*m_biomarkerOccPerChann[channIdx].spindleRippleRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpindleRipple / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spindleFrRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpindleFR / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;
	m_biomarkerOccPerChann[channIdx].spindleHfoRatePerMinute = (double)m_biomarkerOccPerChann[channIdx].totalSpindleHFO / (double)m_biomarkerOccPerChann[channIdx].nrMinutes;*/
}

bool HFO_Detector::generateFileBiomarkersStatistics() {
	
	std::string outputPath = m_outputDirectory + "/MOSSDET_Output/";
	std::string patPath = outputPath + m_patientName;
	std::string filename = patPath + "/" + m_patientName + "_BiomarkerOcurrenceRate.txt";
	std::ofstream biomarkerRatesFile;

	if (!fileExists(filename.c_str())) {
		//Header
		biomarkerRatesFile.open(filename);
		biomarkerRatesFile << "Channel\t" << "DateAndTime\t" << "DayNr\t" << "MinNr\t" << "MinsAvg\t" << "Ripples/Min\t" << "FR/Min\t" << "HFO/Min\t" << "IES/Min\t" << "IES-Ripples/Min\t" << "IES-FR/Min\t" << "IES-HFO/Min\t" << "SS-Ripples/Min\t" << "SS-FR/Min\t" << "SS-HFO/Min\n";
	}
	else {
		biomarkerRatesFile.open(filename, std::ios::app);    // open file for appending
	}

	/*std::string firstDateString;
	struct tm firstDate = { 0, 0, 0 };  // nominal time midday (arbitrary).
	firstDate.tm_year = m_fileStartDateAndTime.year - 1900;
	firstDate.tm_mon = m_fileStartDateAndTime.month - 1;  // note: zero indexed
	firstDate.tm_mday = m_fileStartDateAndTime.day;       // note: not zero indexed
	firstDate.tm_hour = m_fileStartDateAndTime.hours;     // note: zero indexed
	firstDate.tm_min = m_fileStartDateAndTime.minutes;       // note: not zero indexed
	firstDate.tm_sec = m_fileStartDateAndTime.seconds;       // note: not zero indexed

	time_t first_date_seconds = mktime(&firstDate);
	firstDate = *localtime(&first_date_seconds);
	firstDateString = asctime(&firstDate);

	struct tm newDate = { 0, 0, 0 };  // nominal time midday (arbitrary).
	time_t new_date_seconds = mktime(&firstDate) + (m_minutesAnalyzed) * 60;
	newDate = *localtime(&new_date_seconds);*/
	std::string dateAndTimeString = "d";// std::to_string(newDate.tm_year) + "-" + std::to_string(newDate.tm_mon) + "-" + std::to_string(newDate.tm_mday) + "@" + std::to_string(newDate.tm_hour) + ":" + std::to_string(newDate.tm_min) + ":" + std::to_string(newDate.tm_sec);
	int dayNr = ((m_minutesAnalyzed/60.0) / 24);

	for (int ch = 0; ch < m_nrChannels; ++ch) {
		std::string channelName = m_biomarkerOccPerChann[ch].channelName;
		channelName.erase(std::remove(channelName.begin(), channelName.end(), '_'), channelName.end());
		biomarkerRatesFile << channelName << "\t";
		biomarkerRatesFile << dateAndTimeString << "\t";
		biomarkerRatesFile << dayNr << "\t";
		biomarkerRatesFile << m_minutesAnalyzed << "\t";

		biomarkerRatesFile << m_biomarkerOccPerChann[ch].nrMinutes << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].rippleRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].frRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].hfoRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spikeRatePerMinute << "\t";

		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spikeAndRippleRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spikeAndFrRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spikeAndHfoRatePerMinute << "\t";

		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spindleRippleRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spindleFrRatePerMinute << "\t";
		biomarkerRatesFile << m_biomarkerOccPerChann[ch].spindleHfoRatePerMinute << "\n";
		biomarkerRatesFile.flush();
		if (biomarkerRatesFile.bad())
			return false;
	}
	biomarkerRatesFile.close();
	if (biomarkerRatesFile.bad())
		return false;

	return true;
}