#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::min_element, std::max_element, std::find
#include <math.h>       /* sqrt */
#include <numeric>
#include <ctime>
#include <direct.h>
#include <bitset>

#include <dlib/svm.h>

#include "DataTypes.h"
#include "SVM_Detector.h"
#include "SignalProcessing.h"
#include "HFO_EventsHandler.h"

#include <omp.h>

using namespace dlib;
typedef matrix<double, 0, 1> sample_type;
typedef radial_basis_kernel<sample_type> kernel_type;
typedef decision_function<kernel_type> dec_funct_type;
typedef normalized_function<dec_funct_type> funct_type;
typedef probabilistic_decision_function<kernel_type> probabilistic_funct_type;
typedef normalized_function<probabilistic_funct_type> pfunct_type;

SVM_Detector::SVM_Detector(std::string eoiName, std::string patientName, std::string functionName, std::string trainingDataFileName, 
	std::vector<unsigned> selectedFeatures, std::string outputPath, bool saveTextFile) {

	// Classifier parameters
	m_trainingDataFileName = trainingDataFileName;
	m_selectedFeats = selectedFeatures;
	m_gamma = 0;
	m_nu = 0;
	m_probThreshold = 0;
	m_epochLength = 0;
	m_savedFunctionName = functionName;
	m_eoiName = eoiName;

	// Characterized Signal parameters
	m_patientName = patientName;
	m_filename;
	m_inputEvents;
	m_nrEvents = 0;
	m_samplingRate = 0;
	m_signalDuration = 0;
	m_avgMarkDuration = 0;
	m_expertAgreement = 0;

	// Detections Results
	m_nrNoOscills_EOI = 0;
	m_nrShort_EOI = 0;
	m_nrLong_EOI = 0;
	m_nrNoHistoryCompliant_EOI = 0;
	m_avgDetectionDuration = 0;
	m_chSpecifificDetecEpochs;
	m_detectedEOI;

	m_svmPerformance.falseNegatives = 0;
	m_svmPerformance.falsePositives = 0;
	m_svmPerformance.trueNegatives = 0;
	m_svmPerformance.truePositives = 0;

	m_saveTextFile = saveTextFile;
	m_outputPath = outputPath + "/MOSSDET_Output";
	m_patPath = outputPath + "/" + m_patientName;
	m_patDetectionFilesPath = m_outputPath + "/" + m_patientName + "/DetectionFiles";
	_mkdir(m_outputPath.c_str());
	//_mkdir(m_patPath.c_str());
	_mkdir(m_patDetectionFilesPath.c_str());
}

SVM_Detector::~SVM_Detector() {
}

bool SVM_Detector::detectEpochs(double chNr) {

	unsigned numTestingEvs = m_inputEvents.size();

	bool rippleDetection = m_eoiName.find("ast") == std::string::npos;

	funct_type learned_function;

	unsigned nrSVs = 0;
	//std::cout << "deserialize" << m_savedFunctionName << std::endl;
	deserialize(m_savedFunctionName) >> learned_function;
	nrSVs = learned_function.function.basis_vectors.size();
	//std::cout << "decision function read ok" << std::endl;

	// Read samples and learn teh mean and standard deviation in order to provide a normalizer to the SVM
	std::vector<sample_type> samples;
	vector_normalizer<sample_type> normalizer;
	unsigned nrSelectedFeatures = m_selectedFeats.size();
	for (unsigned evIdx = 0; evIdx < numTestingEvs; ++evIdx) {
		sample_type samp;
		samp.set_size(nrSelectedFeatures);
		unsigned numTrainingFeatures = m_inputEvents[evIdx].inputVals.size();
		unsigned nrReadFeatures = 0;
		for (unsigned featIdx = 0; featIdx < numTrainingFeatures; ++featIdx) {
			if (std::find(m_selectedFeats.begin(), m_selectedFeats.end(), featIdx) != m_selectedFeats.end()) {
				double val = m_inputEvents[evIdx].inputVals[featIdx];
				if (isnan(val) || isinf(val)) {
					val = 0;
				}

				samp(nrReadFeatures) = val;
				nrReadFeatures++;
			}
		}
		samples.push_back(samp);
	}
	normalizer.train(samples);		//here the mean and standardr deviation of the provided samples are learned

	m_epochLength = 0;
	m_chSpecifificDetecEpochs[chNr].clear();
	m_chSpecifificDetecEpochs[chNr].resize(numTestingEvs);
	#pragma omp parallel for
	for (long evIdx = 0; evIdx < numTestingEvs; ++evIdx) {
		funct_type svm_binaryClassifier;
		svm_binaryClassifier.normalizer = normalizer;
		svm_binaryClassifier.function = learned_function.function;
		double result = svm_binaryClassifier(samples[evIdx]);
		bool epochVal = result >= m_probThreshold;
		m_chSpecifificDetecEpochs[chNr][evIdx] = (epochVal);
		if (evIdx > 0) {
			m_epochLength += (m_inputEvents[evIdx].time - m_inputEvents[evIdx - 1].time) * 2;
		}
	}
	m_epochLength /= (numTestingEvs-1);

	return true;
}

bool SVM_Detector::getEEG_EventsFromDetectedEpochs(unsigned sharedChs, std::string channName) {

	double detectionOverlap = 0;

	unsigned nrChs = m_chSpecifificDetecEpochs.size();
	unsigned nrEpochs = m_chSpecifificDetecEpochs[0].size();

	unsigned nrShortEOIs = 0;
	unsigned nrLongEOIs = 0;

	double eventMaxAmplitude = 0;
	double eventMaxPower = 0.0;
	double eventMaxSpectralPeak = 0.0;
	std::vector<double> backgroundAmplitude;
	std::vector<double> backgroundPower;

	//Ripple values for min/max duration and feature indices
	double minEOI_Duration = 0.0;// m_epochLength
	double maxEOI_Duration = 0.2;

	int amplitudeFeatIdx = 40;
	int powerFeatIdx = 43;
	int spectrPeakFeatIdx = 6;

	bool fastRipple = m_eoiName.find("ast") != std::string::npos;
	bool spike = m_eoiName.find("ike") != std::string::npos;
	bool spindle = m_eoiName.find("indle") != std::string::npos;

	if (fastRipple) {
		minEOI_Duration = 0.0;// 0.035;
		maxEOI_Duration = 0.2;

		amplitudeFeatIdx = 54;
		powerFeatIdx = 57;
		spectrPeakFeatIdx = 10;
	}

	if (spike) {
		minEOI_Duration = 0.045;// 0.035; 0.06;
		maxEOI_Duration = 0.35;

		amplitudeFeatIdx = 26;
		powerFeatIdx = 29;
		spectrPeakFeatIdx = 2;
	}

	if (spindle) {
		minEOI_Duration = 0.45;// 0.035;
		maxEOI_Duration = 3.5;

		amplitudeFeatIdx = 31;
		powerFeatIdx = 35;
		spectrPeakFeatIdx = 36;		// for spindles this is nr. ZeroCrossings
	}

	bool eoiEpoch = false, eoiStart = false;

	unsigned channel = 0;
	double startTime = 0, endTime = 0, duration = 0;
	long startEvIdx = 0, endEvIdx = 0;

	//writeFeaturesAndOutput(channName);

	for (unsigned evIdx = 1; evIdx < nrEpochs; ++evIdx) {
		eoiEpoch = m_chSpecifificDetecEpochs[0][evIdx];
		// get background feature averages
		if (!eoiEpoch && !eoiStart) {
			backgroundAmplitude.push_back(m_inputEvents[evIdx].inputVals[amplitudeFeatIdx]);
			backgroundPower.push_back(m_inputEvents[evIdx].inputVals[powerFeatIdx]);

			if (backgroundAmplitude.size() *  m_epochLength > 1) {		//refresh vectors when tehy are larger than one minute
				backgroundAmplitude.erase(backgroundAmplitude.begin());
				backgroundPower.erase(backgroundPower.begin());
			}
		}
		// Event Start
		if (eoiEpoch && !eoiStart) {
			startTime = m_inputEvents.at(evIdx).time - m_epochLength;
			eoiStart = true;
			startEvIdx = evIdx;
		}
		// Event has started, characterize it
		if (eoiEpoch && eoiStart) {
			eventMaxAmplitude = m_inputEvents[evIdx].inputVals[amplitudeFeatIdx] > eventMaxAmplitude ? m_inputEvents[evIdx].inputVals[amplitudeFeatIdx] : eventMaxAmplitude;
			eventMaxPower = m_inputEvents[evIdx].inputVals[powerFeatIdx] > eventMaxPower ? m_inputEvents[evIdx].inputVals[powerFeatIdx] : eventMaxPower;
			eventMaxSpectralPeak = m_inputEvents[evIdx].inputVals[spectrPeakFeatIdx] > eventMaxSpectralPeak ? m_inputEvents[evIdx].inputVals[spectrPeakFeatIdx] : eventMaxSpectralPeak;
		}
		// Event end
		if (!eoiEpoch && eoiStart) {
			eoiStart = false;
			channel = m_inputEvents.at(evIdx).channel;
			endTime = m_inputEvents.at(evIdx).time - m_epochLength;
			duration = endTime - startTime;
			endEvIdx = evIdx;

			EOI_Event newDetection;
			newDetection.channelName = channName;
			newDetection.description = m_eoiName;
			newDetection.channelNr = channel;
			newDetection.startTime = startTime;
			newDetection.endTime = endTime;
			newDetection.duration = duration;
			//detection features
			newDetection.amplitude = eventMaxAmplitude;
			newDetection.power = eventMaxPower;
			newDetection.spectralPeak = eventMaxSpectralPeak;
			//Background fatures
			if (!backgroundAmplitude.empty()) {
				newDetection.backgroundAmplitude = std::accumulate(backgroundAmplitude.begin(), backgroundAmplitude.end(), 0.0) / backgroundAmplitude.size();
				newDetection.backgroundPower = std::accumulate(backgroundPower.begin(), backgroundPower.end(), 0.0) / backgroundPower.size();
				SignalProcessing varianceCalc;
				varianceCalc.setVarianceData(backgroundPower);
				double backgroundVar = varianceCalc.runVariance(0, backgroundPower.size());
				newDetection.backgroundStdDev = sqrt(backgroundVar);
			}

			//reset event features
			eventMaxAmplitude = 0;
			eventMaxPower = 0;
			eventMaxSpectralPeak = 0;

			// discard detections because of min or max duration
			if (duration < minEOI_Duration) {
				m_nrShort_EOI++;
				continue;
			}

			if (duration > maxEOI_Duration) {
				//if (duration > 7)
				m_nrLong_EOI++;
				endTime = startTime + maxEOI_Duration;
				duration = maxEOI_Duration;
				continue;
			}

			// save detection
			unsigned nrEOIs = m_detectedEOI.size();
			double minInterDetectDistance = m_epochLength;
			if (spike) {
				minInterDetectDistance = 0.25;
			}
			if (spindle) {
				minInterDetectDistance = 0.45;
			}

			if (nrEOIs > 0) {
				if (newDetection.description.compare(m_detectedEOI[nrEOIs - 1].description) == 0) {
					if (startTime - m_detectedEOI[nrEOIs - 1].endTime < minInterDetectDistance) {//Join detections
						if (endTime - m_detectedEOI[nrEOIs - 1].startTime > maxEOI_Duration) {// Truncate joined detectios which are too long
							m_detectedEOI[nrEOIs - 1].endTime = m_detectedEOI[nrEOIs - 1].startTime + maxEOI_Duration;
						}
						else {
							m_detectedEOI[nrEOIs - 1].endTime = endTime;
						}
						m_detectedEOI[nrEOIs - 1].amplitude = m_detectedEOI[nrEOIs - 1].amplitude > newDetection.amplitude ? m_detectedEOI[nrEOIs - 1].amplitude : newDetection.amplitude; // choose highest amplitude from teh two devents to join
						m_detectedEOI[nrEOIs - 1].power = m_detectedEOI[nrEOIs - 1].power > newDetection.power ? m_detectedEOI[nrEOIs - 1].power : newDetection.power; // choose highest power from teh two devents to join
						m_detectedEOI[nrEOIs - 1].spectralPeak = m_detectedEOI[nrEOIs - 1].spectralPeak > newDetection.spectralPeak ? m_detectedEOI[nrEOIs - 1].spectralPeak : newDetection.spectralPeak; // choose highest spectral peak from teh two devents to join
						m_detectedEOI[nrEOIs - 1].backgroundAmplitude = m_detectedEOI[nrEOIs - 1].backgroundAmplitude;
						m_detectedEOI[nrEOIs - 1].backgroundPower = m_detectedEOI[nrEOIs - 1].backgroundPower;
						m_detectedEOI[nrEOIs - 1].backgroundStdDev = m_detectedEOI[nrEOIs - 1].backgroundStdDev;
					}
					else {
						m_detectedEOI.push_back(newDetection);
					}
				}
				else {
					m_detectedEOI.push_back(newDetection);
				}
			}
			else {
				m_detectedEOI.push_back(newDetection);
			}
		}
	}
	return true;
}

bool SVM_Detector::generateDetectionsAndFeaturesFile(double onsetMinDiff, double minSharedPercent) {

	std::string filename = m_patDetectionFilesPath + "/" + m_patientName + "_" + m_eoiName + "DetectionsAndFeatures.txt";

	double nrDetections = m_detectedEOI.size();
	std::vector<bool> touchedDetects(nrDetections, false);
	m_avgDetectionDuration = 0;

	std::ofstream detectsFile;
	if (!fileExists(filename.c_str())) {
		detectsFile.open(filename);
		//Header
		detectsFile << "Description\t" << "ChannelName\t" << "StartTime(s)\t" << "EndTime(s)\t" << "MaxEventAmplitude\t" << "MaxEventPower\t" << "MaxEventSpectralPeak (Hz)\t" << "AvgBackgroundAmplitude\t" << "AvgBackgroundPower\t" << "BackgroundStdDev\n";
	}
	else {
		detectsFile.open(filename, std::ios::app);    // open file for appending
	}
	detectsFile.precision(16);

	for (long di = 0; di < nrDetections; di++) {
		detectsFile << m_detectedEOI[di].description << "\t" << m_detectedEOI[di].channelName << "\t" << m_detectedEOI[di].startTime << "\t" << m_detectedEOI[di].endTime << "\t"
			<< m_detectedEOI[di].amplitude << "\t" << m_detectedEOI[di].power << "\t" << m_detectedEOI[di].spectralPeak << "\t"
			<< m_detectedEOI[di].backgroundAmplitude << "\t" << m_detectedEOI[di].backgroundPower << "\t" << m_detectedEOI[di].backgroundStdDev << "\n";
		detectsFile.flush();
	}
	detectsFile.close();
	return true;
}

bool SVM_Detector::fileExists(const char *fileName) {
	std::ifstream infile(fileName);
	return infile.good();
}


bool SVM_Detector::getEvents(std::vector<EOI_Event> &readEvents) {
	readEvents = m_detectedEOI;
	return true;
}
