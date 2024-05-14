#pragma once

#include <string>
#include <vector>
#include <array>
#include <algorithm>    // std::min_element, std::max_element
#include <numeric> 

#include "DataTypes.h"
#include "SignalProcessing.h"

#define CONSIDER_PROPAGATION_HFO		false
#define CONSIDER_PROPAGATION_SPIKES		false
#define PROPAGATED_ANNOT_LABEL		0.5
#define RIPPLE_LABEL				1.0
#define FAST_RIPPLE_LABEL			2.0
#define RIPPLE_FR_LABEL				3.0
#define SPIKE_LABEL					10

class HFO_EventsHandler {
public:
	HFO_EventsHandler(){};
};

struct EOI_Features {

	EOI_Features() {
		//Frequency Features Wavelet
		featureNames.push_back("m_sumEpochSubHFOPower");
		featureNames.push_back("m_avgEpochMaxSubHFOPower");
		featureNames.push_back("m_avgEpochSptrlPeakSubHFO");
		featureNames.push_back("m_avgEpochMinMaxRatioSubHFO");

		featureNames.push_back("m_sumEpochRpplPower");
		featureNames.push_back("m_avgEpochMaxRpplPower");
		featureNames.push_back("m_avgEpochSpctrlPeakRipple");
		featureNames.push_back("m_avgEpochMinMaxRatioRipple");

		featureNames.push_back("m_sumEpochFastRpplPower");
		featureNames.push_back("m_avgEpochMaxFastRpplPower");
		featureNames.push_back("m_avgEpochSpectralPeakFR");
		featureNames.push_back("m_avgEpochMinMaxRatioFR");

		// Raw Features
		//Amplitude Features
		featureNames.push_back("m_maxAmplRaw");
		featureNames.push_back("m_amplVarianceRaw");
		featureNames.push_back("m_lineLengthRaw");
		featureNames.push_back("m_powerRaw");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrRaw");
		featureNames.push_back("m_peaksRateNrRaw");
		featureNames.push_back("m_symmetryRaw");
		featureNames.push_back("m_asymmetryRaw");
		// Hjorth Features
		featureNames.push_back("m_mobilityRaw");
		featureNames.push_back("m_complexityRaw");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyRaw");
		//Average Derivative
		featureNames.push_back("m_sumFirstDerivRaw");
		featureNames.push_back("m_sumSecondDerivRaw");
		featureNames.push_back("m_sumThirdDerivRaw");

		//NoDC + Notch filtered features
		//Amplitude Features
		featureNames.push_back("m_maxAmplNotchedDC");
		featureNames.push_back("m_amplVarianceNotchedDC");
		featureNames.push_back("m_lineLengthNotchedDC");
		featureNames.push_back("m_powerNotchedDC");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrNotchedDC");
		featureNames.push_back("m_peaksRateNrNotchedDC");
		featureNames.push_back("m_symmetryNotchedDC");
		featureNames.push_back("m_asymmetryNotchedDC");
		// Hjorth Features
		featureNames.push_back("m_mobilityNotchedDC");
		featureNames.push_back("m_complexityNotchedDC");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyNotchedDC");
		//Average Derivative
		featureNames.push_back("m_sumFirstDerivNotchedDC");
		featureNames.push_back("m_sumSecondDerivNotchedDC");
		featureNames.push_back("m_sumThirdDerivNotchedDC");

		//Ripple Bandpassed Features
		//Amplitude Features
		featureNames.push_back("m_maxAmplRppl");
		featureNames.push_back("m_amplVarianceRppl");
		featureNames.push_back("m_lineLengthRppl");
		featureNames.push_back("m_powerRppl");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrRppl");
		featureNames.push_back("m_peaksRateNrRppl");
		featureNames.push_back("m_symmetryRppl");
		featureNames.push_back("m_asymmetryRppl");
		// Hjorth Features
		featureNames.push_back("m_mobilityRppl");
		featureNames.push_back("m_complexityRppl");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyRppl");
		//Average Derivative
		featureNames.push_back("m_sumFirstDerivRppl");
		featureNames.push_back("m_sumSecondDerivRppl");
		featureNames.push_back("m_sumThirdDerivRppl");

		//Fast-Ripple Bandpassed Features
		//Amplitude Features
		featureNames.push_back("m_maxAmplFastRppl");
		featureNames.push_back("m_amplVarianceFastRppl");
		featureNames.push_back("m_lineLengthFastRppl");
		featureNames.push_back("m_powerFastRppl");
		//Waveform Features
		featureNames.push_back("m_zeroCrossingsNrFastRppl");
		featureNames.push_back("m_peaksRateNrFastRppl");
		featureNames.push_back("m_symmetryFastRppl");
		featureNames.push_back("m_asymmetryFastRppl");
		// Hjorth Features
		featureNames.push_back("m_mobilityFastRppl");
		featureNames.push_back("m_complexityFastRppl");
		//Teager Energy Features	
		featureNames.push_back("m_teagerEnergyFastRppl");
		//Average Derivative
		featureNames.push_back("m_sumFirstDerivFastRppl");
		featureNames.push_back("m_sumSecondDerivFastRppl");
		featureNames.push_back("m_sumThirdDerivFastRppl");

		//Frequency Features Wavelet
		features.push_back(&m_sumEpochSubHFOPower);
		features.push_back(&m_avgEpochMaxSubHFOPower);
		features.push_back(&m_avgEpochSptrlPeakSubHFO);
		features.push_back(&m_avgEpochMinMaxRatioSubHFO);

		features.push_back(&m_sumEpochRpplPower);
		features.push_back(&m_avgEpochMaxRpplPower);
		features.push_back(&m_avgEpochSpctrlPeakRipple);
		features.push_back(&m_avgEpochMinMaxRatioRipple);

		features.push_back(&m_sumEpochFastRpplPower);
		features.push_back(&m_avgEpochMaxFastRpplPower);
		features.push_back(&m_avgEpochSpectralPeakFR);
		features.push_back(&m_avgEpochMinMaxRatioFR);


		// Raw Features
		//Amplitude Features
		features.push_back(&m_maxAmplRaw);
		features.push_back(&m_amplVarianceRaw);
		features.push_back(&m_lineLengthRaw);
		features.push_back(&m_powerRaw);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrRaw);
		features.push_back(&m_peaksRateNrRaw);
		features.push_back(&m_symmetryRaw);
		features.push_back(&m_asymmetryRaw);
		// Hjorth Features
		features.push_back(&m_mobilityRaw);
		features.push_back(&m_complexityRaw);
		//Teager Energy Features	
		features.push_back(&m_teagerEnergyRaw);
		//Average Derivative
		features.push_back(&m_sumFirstDerivRaw);
		features.push_back(&m_sumSecondDerivRaw);
		features.push_back(&m_sumThirdDerivRaw);

		//NoDC + Notch filtered features
		//Amplitude Features
		features.push_back(&m_maxAmplNotchedDC);
		features.push_back(&m_amplVarianceNotchedDC);
		features.push_back(&m_lineLengthNotchedDC);
		features.push_back(&m_powerNotchedDC);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrNotchedDC);
		features.push_back(&m_peaksRateNrNotchedDC);
		features.push_back(&m_symmetryNotchedDC);
		features.push_back(&m_asymmetryNotchedDC);
		// Hjorth Features
		features.push_back(&m_mobilityNotchedDC);
		features.push_back(&m_complexityNotchedDC);
		//Teager Energy Features	
		features.push_back(&m_teagerEnergyNotchedDC);
		//Average Derivative
		features.push_back(&m_sumFirstDerivNotchedDC);
		features.push_back(&m_sumSecondDerivNotchedDC);
		features.push_back(&m_sumThirdDerivNotchedDC);

		//Ripple Bandpassed Features
		//Amplitude Features
		features.push_back(&m_maxAmplRppl);
		features.push_back(&m_amplVarianceRppl);
		features.push_back(&m_lineLengthRppl);
		features.push_back(&m_powerRppl);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrRppl);
		features.push_back(&m_peaksRateNrRppl);
		features.push_back(&m_symmetryRppl);
		features.push_back(&m_asymmetryRppl);
		// Hjorth Features
		features.push_back(&m_mobilityRppl);
		features.push_back(&m_complexityRppl);
		// Other Feats
		features.push_back(&m_teagerEnergyRppl);
		//Average Derivative
		features.push_back(&m_sumFirstDerivRppl);
		features.push_back(&m_sumSecondDerivRppl);
		features.push_back(&m_sumThirdDerivRppl);

		//Fast-Ripple Bandpassed Features
		//Amplitude Features
		features.push_back(&m_maxAmplFastRppl);
		features.push_back(&m_amplVarianceFastRppl);
		features.push_back(&m_lineLengthFastRppl);
		features.push_back(&m_powerFastRppl);
		//Waveform Features
		features.push_back(&m_zeroCrossingsNrFastRppl);
		features.push_back(&m_peaksRateNrFastRppl);
		features.push_back(&m_symmetryFastRppl);
		features.push_back(&m_asymmetryFastRppl);
		// Hjorth Features
		features.push_back(&m_mobilityFastRppl);
		features.push_back(&m_complexityFastRppl);
		// Other Feats
		features.push_back(&m_teagerEnergyFastRppl);
		//Average Derivative
		features.push_back(&m_sumFirstDerivFastRppl);
		features.push_back(&m_sumSecondDerivFastRppl);
		features.push_back(&m_sumThirdDerivFastRppl);
	}

	std::vector<std::vector<std::vector<double>>*> features;			//feature, channel, epoch
	std::vector<std::string> featureNames;

	//Frequency Features Wavelet[0:8]
	std::vector<std::vector<double>> m_sumEpochSubHFOPower;				//0
	std::vector<std::vector<double>> m_avgEpochMaxSubHFOPower;			//1
	std::vector<std::vector<double>> m_avgEpochSptrlPeakSubHFO;			//2
	std::vector<std::vector<double>> m_avgEpochMinMaxRatioSubHFO;		//3

	std::vector<std::vector<double>> m_sumEpochRpplPower;				//4
	std::vector<std::vector<double>> m_avgEpochMaxRpplPower;			//5
	std::vector<std::vector<double>> m_avgEpochSpctrlPeakRipple;		//6
	std::vector<std::vector<double>> m_avgEpochMinMaxRatioRipple;		//7

	std::vector<std::vector<double>> m_sumEpochFastRpplPower;			//8
	std::vector<std::vector<double>> m_avgEpochMaxFastRpplPower;		//9
	std::vector<std::vector<double>> m_avgEpochSpectralPeakFR;			//10
	std::vector<std::vector<double>> m_avgEpochMinMaxRatioFR;			//11

	// Raw Signal Features[9:21]
	//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplRaw;						//12
	std::vector<std::vector<double>> m_amplVarianceRaw;					//13
	std::vector<std::vector<double>> m_lineLengthRaw;					//14
	std::vector<std::vector<double>> m_powerRaw;						//15
																		//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrRaw;				//16
	std::vector<std::vector<double>> m_peaksRateNrRaw;					//17
	std::vector<std::vector<double>> m_symmetryRaw;						//18
	std::vector<std::vector<double>> m_asymmetryRaw;					//19
																		// Hjorth Features
	std::vector<std::vector<double>> m_mobilityRaw;						//20
	std::vector<std::vector<double>> m_complexityRaw;					//21
																		//Teager Energy Features	
	std::vector<std::vector<double>> m_teagerEnergyRaw;					//22
																		//Average Derivative	
	std::vector<std::vector<double>> m_sumFirstDerivRaw;				//23
	std::vector<std::vector<double>> m_sumSecondDerivRaw;				//24
	std::vector<std::vector<double>> m_sumThirdDerivRaw;				//25

																		//NoDC + Notch filtered signal features[22:32]
																		//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplNotchedDC;				//26
	std::vector<std::vector<double>> m_amplVarianceNotchedDC;			//27
	std::vector<std::vector<double>> m_lineLengthNotchedDC;				//28
	std::vector<std::vector<double>> m_powerNotchedDC;					//29
																		//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrNotchedDC;		//30
	std::vector<std::vector<double>> m_peaksRateNrNotchedDC;			//31
	std::vector<std::vector<double>> m_symmetryNotchedDC;				//32
	std::vector<std::vector<double>> m_asymmetryNotchedDC;				//33
																		// Hjorth Features
	std::vector<std::vector<double>> m_mobilityNotchedDC;				//34
	std::vector<std::vector<double>> m_complexityNotchedDC;				//35
																		//Teager Energy Features
	std::vector<std::vector<double>> m_teagerEnergyNotchedDC;			//36
																		//Average Derivative	
	std::vector<std::vector<double>> m_sumFirstDerivNotchedDC;			//37
	std::vector<std::vector<double>> m_sumSecondDerivNotchedDC;			//38
	std::vector<std::vector<double>> m_sumThirdDerivNotchedDC;			//39

																		//Ripple Bandpassed Features[33:43]
																		//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplRppl;						//40
	std::vector<std::vector<double>> m_amplVarianceRppl;				//41
	std::vector<std::vector<double>> m_lineLengthRppl;					//42
	std::vector<std::vector<double>> m_powerRppl;						//43
																		//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrRppl;				//44
	std::vector<std::vector<double>> m_peaksRateNrRppl;					//45
	std::vector<std::vector<double>> m_symmetryRppl;					//46
	std::vector<std::vector<double>> m_asymmetryRppl;					//47
																		// Hjorth Features
	std::vector<std::vector<double>> m_mobilityRppl;					//48
	std::vector<std::vector<double>> m_complexityRppl;					//49
																		//Teager Energy Features	
	std::vector<std::vector<double>> m_teagerEnergyRppl;				//50
																		//Average Derivative	
	std::vector<std::vector<double>> m_sumFirstDerivRppl;				//51
	std::vector<std::vector<double>> m_sumSecondDerivRppl;				//52
	std::vector<std::vector<double>> m_sumThirdDerivRppl;				//53

																		//Fast Ripple Bandpassed Features[44:54]
																		//Amplitude Features
	std::vector<std::vector<double>> m_maxAmplFastRppl;					//54
	std::vector<std::vector<double>> m_amplVarianceFastRppl;			//55
	std::vector<std::vector<double>> m_lineLengthFastRppl;				//56
	std::vector<std::vector<double>> m_powerFastRppl;					//57

																		//Waveform Features
	std::vector<std::vector<double>> m_zeroCrossingsNrFastRppl;			//58
	std::vector<std::vector<double>> m_peaksRateNrFastRppl;				//59
	std::vector<std::vector<double>> m_symmetryFastRppl;				//60
	std::vector<std::vector<double>> m_asymmetryFastRppl;				//61
																		// Hjorth Features
	std::vector<std::vector<double>> m_mobilityFastRppl;				//62
	std::vector<std::vector<double>> m_complexityFastRppl;				//63
																		//Teager Energy Features	
	std::vector<std::vector<double>> m_teagerEnergyFastRppl;			//64
																		//Average Derivative	
	std::vector<std::vector<double>> m_sumFirstDerivFastRppl;			//65
	std::vector<std::vector<double>> m_sumSecondDerivFastRppl;			//66
	std::vector<std::vector<double>> m_sumThirdDerivFastRppl;			//67

	std::vector<double> m_time;
	std::vector<std::vector<double>> m_EOI_Mark;

	void setMatrices(unsigned nrChannels) {
		for (unsigned i = 0; i < features.size(); ++i) {
			features[i]->resize(nrChannels);
		}
		m_EOI_Mark.resize(nrChannels);
	}

	void transformFeaturesSmoothen() {
		for (unsigned featIdx = 0; featIdx < features.size(); ++featIdx) { // Iterate through all features
			unsigned nrChannels = features[featIdx]->size();
			for (unsigned chIdx = 0; chIdx < nrChannels; ++chIdx) {	// Iterate through all channels of a feature
				unsigned nrSamples = features[featIdx]->at(chIdx).size();
				if (nrSamples != 0) {
					unsigned avgSize = 4;
					for (unsigned sample = 0; sample < nrSamples; ++sample) {	// Iterate and normalize all channel samples
						double avg = 0;
						if (sample >= 1 && sample < nrSamples - 1) {
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
