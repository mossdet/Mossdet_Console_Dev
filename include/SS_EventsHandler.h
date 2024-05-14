#pragma once

#include <string>
#include <vector>
#include <array>
#include <algorithm>    // std::min_element, std::max_element
#include <numeric> 

#include "DataTypes.h"
#include "SignalProcessing.h"

class SS_EventsHandler {
public:
	SS_EventsHandler() {};
};
struct SS_Features {

	SS_Features() {
		//Frequency Features Wavelet
		featureNames.push_back("m_subSpindle_PowWavelet");
		featureNames.push_back("m_spindle_PowWavelet");
		featureNames.push_back("m_supraSpindle_PowWavelet");

		// Raw Features
		//Amplitude Features
		featureNames.push_back("m_maxAmplRaw");
		featureNames.push_back("m_amplVarianceRaw");
		featureNames.push_back("m_meanAmplitudeRaw");
		featureNames.push_back("m_lineLengthRaw");
		featureNames.push_back("m_powerRaw");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrRaw");
		featureNames.push_back("m_peaksRateNrRaw");
		featureNames.push_back("m_autocorrelationRaw");
		featureNames.push_back("m_symmetryRaw");
		featureNames.push_back("m_asymmetryRaw");
		// Hjorth Features
		featureNames.push_back("m_mobilityRaw");
		featureNames.push_back("m_complexityRaw");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyRaw");
		featureNames.push_back("m_teagerEnergyAutocorrRaw");

		//NoDC + Notch filtered features
		//Amplitude Features
		featureNames.push_back("m_maxAmplNotchedDC");
		featureNames.push_back("m_amplVarianceNotchedDC");
		featureNames.push_back("m_meanAmplitudeNotcheDC");
		featureNames.push_back("m_lineLengthNotchedDC");
		featureNames.push_back("m_powerNotchedDC");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrNotchedDC");
		featureNames.push_back("m_peaksRateNrNotchedDC");
		featureNames.push_back("m_autocorrelationNotchedDC");
		featureNames.push_back("m_symmetryNotchedDC");
		featureNames.push_back("m_asymmetryNotchedDC");
		// Hjorth Features
		featureNames.push_back("m_mobilityNotchedDC");
		featureNames.push_back("m_complexityNotchedDC");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyNotchedDC");
		featureNames.push_back("m_teagerEnergyAutocorrNotchedDC");

		//Bandpassed Features
		//Amplitude Features
		featureNames.push_back("m_maxAmplBP");
		featureNames.push_back("m_amplVarianceBP");
		featureNames.push_back("m_meanAmplitudeBP");
		featureNames.push_back("m_lineLengthBP");
		featureNames.push_back("m_powerBP");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrBP");
		featureNames.push_back("m_peaksRateNrBP");
		featureNames.push_back("m_autocorrelationBP");
		featureNames.push_back("m_symmetryBP");
		featureNames.push_back("m_asymmetryBP");
		// Hjorth Features
		featureNames.push_back("m_mobilityBP");
		featureNames.push_back("m_complexityBP");
		// Other Feats
		featureNames.push_back("m_teagerEnergyBP");
		featureNames.push_back("m_teagerEnergyAutocorrBP");

		//Frequency Features Wavelet
		features.push_back(&m_subSpindle_PowWavelet);
		features.push_back(&m_spindle_PowWavelet);
		features.push_back(&m_supraSpindle_PowWavelet);

		// Raw Features
		//Amplitude Features
		features.push_back(&m_maxAmplRaw);
		features.push_back(&m_amplVarianceRaw);
		features.push_back(&m_meanAmplitudeRaw);
		features.push_back(&m_lineLengthRaw);
		features.push_back(&m_powerRaw);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrRaw);
		features.push_back(&m_peaksRateNrRaw);
		features.push_back(&m_autocorrelationRaw);
		features.push_back(&m_symmetryRaw);
		features.push_back(&m_asymmetryRaw);
		// Hjorth Features
		features.push_back(&m_mobilityRaw);
		features.push_back(&m_complexityRaw);
		//Teager Energy Features	
		features.push_back(&m_teagerEnergyRaw);
		features.push_back(&m_teagerEnergyAutocorrRaw);

		//NoDC + Notch filtered features
		//Amplitude Features
		features.push_back(&m_maxAmplNotchedDC);
		features.push_back(&m_amplVarianceNotchedDC);
		features.push_back(&m_meanAmplitudeNotchedDC);
		features.push_back(&m_lineLengthNotchedDC);
		features.push_back(&m_powerNotchedDC);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrNotchedDC);
		features.push_back(&m_peaksRateNrNotchedDC);
		features.push_back(&m_autocorrelationNotchedDC);
		features.push_back(&m_symmetryNotchedDC);
		features.push_back(&m_asymmetryNotchedDC);
		// Hjorth Features
		features.push_back(&m_mobilityNotchedDC);
		features.push_back(&m_complexityNotchedDC);
		//Teager Energy Features	
		features.push_back(&m_teagerEnergyNotchedDC);
		features.push_back(&m_teagerEnergyAutocorrNotchedDC);

		//Bandpassed Features
		//Amplitude Features
		features.push_back(&m_maxAmplBP);
		features.push_back(&m_amplVarianceBP);
		features.push_back(&m_meanAmplitudeBP);
		features.push_back(&m_lineLengthBP);
		features.push_back(&m_powerBP);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrBP);
		features.push_back(&m_peaksRateNrBP);
		features.push_back(&m_autocorrelationBP);
		features.push_back(&m_symmetryBP);
		features.push_back(&m_asymmetryBP);
		// Hjorth Features
		features.push_back(&m_mobilityBP);
		features.push_back(&m_complexityBP);
		// Other Feats
		features.push_back(&m_teagerEnergyBP);
		features.push_back(&m_teagerEnergyAutocorrBP);
	}

	std::vector<std::vector<std::vector<double>>*> features;
	std::vector<std::string> featureNames;



	//Frequency Features Wavelet[0:2]
	std::vector<std::vector<double>> m_subSpindle_PowWavelet;			//0
	std::vector<std::vector<double>> m_spindle_PowWavelet;				//1
	std::vector<std::vector<double>> m_supraSpindle_PowWavelet;			//2


																		// Raw Signal Features[3:16]
																		//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplRaw;					//3
	std::vector<std::vector<double>> m_amplVarianceRaw;				//4
	std::vector<std::vector<double>> m_meanAmplitudeRaw;			//5
	std::vector<std::vector<double>> m_lineLengthRaw;				//6
	std::vector<std::vector<double>> m_powerRaw;					//7

																	//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrRaw;			//8
	std::vector<std::vector<double>> m_peaksRateNrRaw;				//9
	std::vector<std::vector<double>> m_autocorrelationRaw;			//10
	std::vector<std::vector<double>> m_symmetryRaw;					//11
	std::vector<std::vector<double>> m_asymmetryRaw;				//12
																	// Hjorth Features
	std::vector<std::vector<double>> m_mobilityRaw;					//13
	std::vector<std::vector<double>> m_complexityRaw;				//14
																	//Teager Energy Features	
	std::vector<std::vector<double>> m_teagerEnergyRaw;				//15
	std::vector<std::vector<double>> m_teagerEnergyAutocorrRaw;		//16

																	//NoDC + Notch filtered signal features[17:30]
																	//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplNotchedDC;			//17
	std::vector<std::vector<double>> m_amplVarianceNotchedDC;		//18
	std::vector<std::vector<double>> m_meanAmplitudeNotchedDC;		//19
	std::vector<std::vector<double>> m_lineLengthNotchedDC;			//20
	std::vector<std::vector<double>> m_powerNotchedDC;				//21
																	//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrNotchedDC;	//22
	std::vector<std::vector<double>> m_peaksRateNrNotchedDC;		//23
	std::vector<std::vector<double>> m_autocorrelationNotchedDC;	//24
	std::vector<std::vector<double>> m_symmetryNotchedDC;			//25
	std::vector<std::vector<double>> m_asymmetryNotchedDC;			//26
																	// Hjorth Features
	std::vector<std::vector<double>> m_mobilityNotchedDC;			//27
	std::vector<std::vector<double>> m_complexityNotchedDC;			//28
																	//Teager Energy Features	
	std::vector<std::vector<double>> m_teagerEnergyNotchedDC;		//29
	std::vector<std::vector<double>> m_teagerEnergyAutocorrNotchedDC;//30

																	 //Bandpassed Features[31:44]
																	 //Amplitude Features
	std::vector<std::vector<double>> m_maxAmplBP;				//31
	std::vector<std::vector<double>> m_amplVarianceBP;			//32
	std::vector<std::vector<double>> m_meanAmplitudeBP;			//33
	std::vector<std::vector<double>> m_lineLengthBP;			//34
	std::vector<std::vector<double>> m_powerBP;					//35
																//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrBP;		//36
	std::vector<std::vector<double>> m_peaksRateNrBP;			//37
	std::vector<std::vector<double>> m_autocorrelationBP;		//38
	std::vector<std::vector<double>> m_symmetryBP;				//39
	std::vector<std::vector<double>> m_asymmetryBP;				//40
																// Hjorth Features
	std::vector<std::vector<double>> m_mobilityBP;				//41
	std::vector<std::vector<double>> m_complexityBP;			//42
																// Teager Feats
	std::vector<std::vector<double>> m_teagerEnergyBP;			//43
	std::vector<std::vector<double>> m_teagerEnergyAutocorrBP;	//44

	std::vector<double> m_time;
	std::vector<std::vector<double>> m_SS_Mark;

	void setMatrices(unsigned nrChannels) {
		for (unsigned i = 0; i < features.size(); ++i) {
			features[i]->resize(nrChannels);
		}
		m_SS_Mark.resize(nrChannels);
	}

	void transformFeaturesSmoothen() {
		for (unsigned featIdx = 0; featIdx < features.size(); ++featIdx) { // Iterate through all features
			size_t nrChannels = features[featIdx]->size();
			for (unsigned chIdx = 0; chIdx < nrChannels; ++chIdx) {	// Iterate through all channels of a feature
				size_t nrSamples = features[featIdx]->at(chIdx).size();
				if (nrSamples != 0) {
					unsigned avgSize = 4;
					for (unsigned sample = 0; sample < nrSamples; ++sample) {	// Iterate and normalize all channel samples
						double avg = 0;
						if (sample > 0 && sample < nrSamples - 1) {
							avg =
								features[featIdx]->at(chIdx).at(sample - 1) +
								features[featIdx]->at(chIdx).at(sample) +
								features[featIdx]->at(chIdx).at(sample + 1);
							avg = avg / 3;
							features[featIdx]->at(chIdx).at(sample) = avg;
						}
					}
				}
			}
		}
	}

	void getMax(unsigned featIdx, unsigned channIdx) {
		if (features[featIdx]->at(channIdx).size() != 0) {
			double max = *std::max_element(features[featIdx]->at(channIdx).begin(), features[featIdx]->at(channIdx).end());
		}
	}

	void normalizeFeatures() {
		for (unsigned i = 0; i < features.size(); ++i) { // Iterate through all features
			for (unsigned chann = 0; chann < features[i]->size(); ++chann) {	// Iterate through all channels of a feature
				if (features[i]->at(chann).size() != 0) {
					vectorStats featTracseStats;
					featTracseStats = SignalProcessing::vectorStatistics(features[i]->at(chann));
					for (unsigned sample = 0; sample < features[i]->at(chann).size(); ++sample) {	// Iterate and normalize all channel samples
						double val = features[i]->at(chann).at(sample);
						features[i]->at(chann).at(sample) = (val - featTracseStats.mean) / featTracseStats.stdrdDev;
					}
				}
			}
		}
	}
};