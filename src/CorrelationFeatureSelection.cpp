#include <fstream>
#include <direct.h>
#include <time.h>

#include "HFO_Detector.h"
#include "SignalProcessing.h"
#include "HFO_EventsHandler.h"
#include "CorrelationFeatureSelection.h"
#include "ancorrelation.h"

CorrelationFeatureSelection::CorrelationFeatureSelection(unsigned nrChannels, std::string patientPath, std::string patientName) {
	m_nrChannels = nrChannels;
	m_multichannel = false;
	m_patientName = patientName;
	m_patientPath = patientPath;
}



bool CorrelationFeatureSelection::getEOI_MarkFeatureCorrelation(bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes, EOI_Features &epochFeatures) {

	/* initialize random seed: */
	srand(time(NULL));

	std::ofstream trainingFile;

	std::string pathName = m_patientPath + "/TrainingFiles/CFS/";
	_mkdir(pathName.c_str());

	std::string filename = m_patientPath + "/TrainingFiles/CFS/" + m_patientName;

	if (m_multichannel)
		filename += "_MCH_Features_Correlation";
	else
		filename += "_SCH_Features_Correlation";

	if (getCFS_Ripple && !getCFS_FastRipple)
		filename += "_Ripple.txt";
	else if (getCFS_FastRipple && !getCFS_Ripple)
		filename += "_FastRipple.txt";
	else if (getCFS_AllHFO)
		filename += "_AllHFO.txt";
	else if (getCFS_SpikeHFO)
		filename += "_SpikeHFO.txt";
	else if (getCFS_SpikesAlone)
		filename += "_SpikesAlone.txt";
	else if (getCFS_AllSpikes)
		filename += "_AllSpike.txt";
	else
		filename += "_error.txt";


	remove(filename.c_str());
	trainingFile.open(filename);
	trainingFile.precision(32);
	trainingFile << "Patient\t" << "Channelnr\t" << "ChannelName\t" << "FeatureNr\t" << "FeatureName\t" << "SecondSignal\t" << "CorrCoeff\n";


	bool validData;
	if (getCFS_Ripple || getCFS_FastRipple || getCFS_AllHFO) {
		validData = relabelDataHFO(getCFS_Ripple, getCFS_FastRipple, getCFS_AllHFO, epochFeatures);
	}
	else if (getCFS_SpikeHFO) {
		validData = relabelDataSpikesHFO(epochFeatures);
	}
	else if (getCFS_SpikesAlone) {
		validData = relabelDataSpikesAlone(epochFeatures);
	}
	else if (getCFS_AllSpikes) {
		validData = relabelDataAllSpikes(epochFeatures);
	}

	bool vectorizationOK = false;
	EOI_Features vectorizedFeatureChannels; 

	if (validData) {
		vectorizationOK = vectorizeChannels(epochFeatures, vectorizedFeatureChannels);
	}

	if(vectorizationOK) {
		m_nrChannels = 1;
		unsigned nrEpochs = vectorizedFeatureChannels.m_time.size();
		unsigned nrFeatures = vectorizedFeatureChannels.features.size();

		//Mark-feature correlation
		std::vector<double> featMarkCorrChannAvg(nrFeatures, 0);
		for (int channelIdx = 0; channelIdx < m_nrChannels; ++channelIdx) {
#pragma omp parallel for
			for (int featIdx = 0; featIdx < nrFeatures; ++featIdx) {
				double corrCoeff;
				unsigned featsSize = vectorizedFeatureChannels.features[featIdx]->at(channelIdx).size();
				if (std::find(m_featsToExclude.begin(), m_featsToExclude.end(), featIdx) == m_featsToExclude.end()) {
					getMaxCrossCorrelationPearson(vectorizedFeatureChannels.features[featIdx]->at(channelIdx), vectorizedFeatureChannels.m_EOI_Mark[channelIdx], corrCoeff, 0, 0);
					//getMaxCrossCorrelationRho(vectorizedFeatureChannels.features[featIdx]->at(channelIdx), vectorizedFeatureChannels.m_EOI_Mark[channelIdx], 0, nrEpochs, corrCoeff, 0, 0);
				}
				else {
					corrCoeff = 0;
				}
				if (!std::isfinite(corrCoeff) || !std::isfinite(corrCoeff) || !std::isfinite(corrCoeff))
					corrCoeff = 0;
				featMarkCorrChannAvg[featIdx] += corrCoeff;
			}
		}
		
#pragma omp parallel for
		for (int featIdx = 0; featIdx < nrFeatures; ++featIdx) {
			featMarkCorrChannAvg[featIdx] /= m_nrChannels;
		}

		//Feature-feature correlation
		std::vector<std::vector<double>> featFeatCorrChannAvg;
		for (unsigned i = 0; i < nrFeatures; ++i) {
			std::vector<double> temp(nrFeatures, 0);
			featFeatCorrChannAvg.push_back(temp);
		}

		for (int channelIdx = 0; channelIdx < m_nrChannels; ++channelIdx) {
			for (int featIdx = 0; featIdx < nrFeatures; ++featIdx) {
#pragma omp parallel for
				for (int compFeatIdx = 0; compFeatIdx < nrFeatures; ++compFeatIdx) {
					double corrCoeff;
					if (std::find(m_featsToExclude.begin(), m_featsToExclude.end(), featIdx) == m_featsToExclude.end()) {
						getMaxCrossCorrelationPearson(vectorizedFeatureChannels.features[featIdx]->at(channelIdx), vectorizedFeatureChannels.features[compFeatIdx]->at(channelIdx), corrCoeff, 0, 0);
						//getMaxCrossCorrelationRho(vectorizedFeatureChannels.features[featIdx]->at(channelIdx), vectorizedFeatureChannels.features[compFeatIdx]->at(channelIdx), 0, nrEpochs, corrCoeff, 0, 0);
					}
					else {
						corrCoeff = 0;
					}
					if (!std::isfinite(corrCoeff) || !std::isfinite(corrCoeff) || !std::isfinite(corrCoeff))
						corrCoeff = 0;
					featFeatCorrChannAvg[featIdx][compFeatIdx] += corrCoeff;
				}
			}
		}

		for (int featIdx = 0; featIdx < nrFeatures; ++featIdx) {
#pragma omp parallel for
			for (int compFeatIdx = 0; compFeatIdx < nrFeatures; ++compFeatIdx) {
				featFeatCorrChannAvg[featIdx][compFeatIdx] /= m_nrChannels;
			}
		}

		getCFS(featMarkCorrChannAvg, featFeatCorrChannAvg, getCFS_Ripple, getCFS_FastRipple, getCFS_AllHFO, getCFS_SpikeHFO, getCFS_SpikesAlone, getCFS_AllSpikes);

		bool writeCorrFile = true;
		if (writeCorrFile) {
			//Write Mark-feature correlation
			for (unsigned featIdx = 0; featIdx < nrFeatures; ++featIdx) {
				double channAvgdCorrCoeff = featMarkCorrChannAvg[featIdx];
				trainingFile << m_patientName << "\t" << m_nrChannels << "\t" << "AllChannels" << "\t" << featIdx << "\t" << vectorizedFeatureChannels.featureNames[featIdx] << "\t" << "EOI_Mark" << "\t" << channAvgdCorrCoeff << "\n";
			}
			//Write Feature-feature correlation
			for (unsigned featIdx = 0; featIdx < nrFeatures; ++featIdx) {
				for (unsigned compFeatIdx = 0; compFeatIdx < nrFeatures; ++compFeatIdx) {
					double channAvgdCorrCoeff = featFeatCorrChannAvg[featIdx][compFeatIdx];
					trainingFile << m_patientName << "\t" << m_nrChannels << "\t" << "AllChannels" << "\t" << featIdx << "\t" << vectorizedFeatureChannels.featureNames[featIdx] << "\t" << vectorizedFeatureChannels.featureNames[compFeatIdx] << "\t" << channAvgdCorrCoeff << "\n";
				}
			}
		}
		
	}
	/*else {
		while (1);
	}*/


	trainingFile.close();

	return true;
}

bool CorrelationFeatureSelection::getCFS(std::vector<double> featMarkCorrChannAvg, std::vector<std::vector<double>> featFeatCorrChannAvg, bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, bool getCFS_SpikeHFO, bool getCFS_SpikesAlone, bool getCFS_AllSpikes) {

	unsigned totNrFeats = featMarkCorrChannAvg.size();
	double featSetMerit = 0;
	double prevFeatSetMerit = 0;
	int prevErasedFeat = 0;
	int prevSize = 0;
	std::vector<unsigned> selectedFeats;
	std::vector<unsigned> posFeatures;
	for (unsigned featIdx = 0; featIdx < totNrFeats; ++featIdx) {
		if (std::find(m_featsToExclude.begin(), m_featsToExclude.end(), featIdx) == m_featsToExclude.end()) {
			selectedFeats.push_back(featIdx);
		}
	}

	unsigned unchangedCount = 0;
	while (unchangedCount < 1000 && selectedFeats.size() > 1) {

		double kFeats = 0;
		double setAvg_MFC = 0;																								//mark-feature correlation averaged across all features of the current subset
		double setAvg_FFC = 0;																								//feature-feature correlation averaged across all features of the current subset

		for (unsigned selFeatIdx = 0; selFeatIdx < selectedFeats.size(); ++selFeatIdx) {

			//Calculate values for merit equation
			kFeats = selectedFeats.size();	

			if (kFeats == 1) {
				break;
			}
																															//Feature-Mark correlation
			for (unsigned traverseIdx = 0; traverseIdx < kFeats; ++traverseIdx) {
				unsigned featIdx = selectedFeats[traverseIdx];
				setAvg_MFC += std::abs(featMarkCorrChannAvg[featIdx]);
			}
			setAvg_MFC /= kFeats;


			//Feature-Feature correlation
			for (unsigned traverseIdxOne = 0; traverseIdxOne < kFeats; ++traverseIdxOne) {
				unsigned featIdxOne = selectedFeats[traverseIdxOne];
				for (unsigned traverseIdxTwo = 0; traverseIdxTwo < kFeats; ++traverseIdxTwo) {
					unsigned featIdxTwo = selectedFeats[traverseIdxTwo];
					setAvg_FFC += std::abs(featFeatCorrChannAvg[featIdxOne][featIdxTwo]);
				}
			}
			setAvg_FFC /= (kFeats*kFeats);
		}

		featSetMerit = (kFeats * setAvg_MFC) / sqrt(kFeats + (kFeats * (kFeats - 1) * setAvg_FFC));
		if (prevFeatSetMerit == 0)
			prevFeatSetMerit = featSetMerit;

		if ((int)(featSetMerit * 1000000) >= (int)(prevFeatSetMerit * 1000000)) {
			int featEraseIdx = (rand() % (int)kFeats);
			int featToErase = selectedFeats[featEraseIdx];

			unsigned searchCount = 0;
			while ((std::find(posFeatures.begin(), posFeatures.end(), featToErase) != posFeatures.end()) && searchCount < 1000) {
				featEraseIdx = (rand() % (int)kFeats);
				featToErase = selectedFeats[featEraseIdx];
				searchCount++;
			}

			if (std::find(posFeatures.begin(), posFeatures.end(), featToErase) == posFeatures.end()) {
				prevErasedFeat = featToErase;
				selectedFeats.erase(selectedFeats.begin() + featEraseIdx);
				prevFeatSetMerit = featSetMerit;
			}
			else
				bool stop = true;

		}
		else if ((int)(featSetMerit * 1000000) < (int)(prevFeatSetMerit * 1000000)) {
			selectedFeats.push_back(prevErasedFeat);
			posFeatures.push_back(prevErasedFeat);
		}

		std::sort(selectedFeats.begin(), selectedFeats.end());
		std::sort(posFeatures.begin(), posFeatures.end());

		if (prevSize != selectedFeats.size()) {
			prevSize = selectedFeats.size();
			unchangedCount = 0;
		}
		else {
			unchangedCount++;
		}
	}

	std::ofstream cfsFile;
	std::string filename = m_patientPath + "/TrainingFiles/CFS/" + m_patientName;
	if (m_multichannel)
		filename += "_CFS_MCH";
	else
		filename += "_CFS_SCH";

	if (getCFS_Ripple)
		filename += "_Ripple.txt";
	else if (getCFS_FastRipple)
		filename += "_FastRipple.txt";
	else if (getCFS_AllHFO)
		filename += "_AllHFO.txt";
	else if (getCFS_SpikeHFO)
		filename += "_SpikeHFO.txt";
	else if (getCFS_SpikesAlone)
		filename += "_SpikesAlone.txt";
	else if (getCFS_AllSpikes)
		filename += "_AllSpikes.txt";
	else
		filename += "_error.txt";

	remove(filename.c_str());
	cfsFile.open(filename);

	EOI_Features feats;
	cfsFile << "SelectedFeatures\n";
	for (unsigned selFeatIdx = 0; selFeatIdx < selectedFeats.size(); ++selFeatIdx) {
		cfsFile << "d:\t" << selectedFeats[selFeatIdx] << "\t" << feats.featureNames[selectedFeats[selFeatIdx]] << "\n";
	}

	cfsFile << "\n\n";
	cfsFile << "PositiveFeatures\n";
	for (unsigned selFeatIdx = 0; selFeatIdx < posFeatures.size(); ++selFeatIdx) {
		cfsFile << posFeatures[selFeatIdx] << "\t" << feats.featureNames[posFeatures[selFeatIdx]] << "\n";
	}
	cfsFile.close();

	return true;
}

long CorrelationFeatureSelection::getMaxCrossCorrelationPearson(std::vector<double> &signalA, std::vector<double> &signalB, double& corrCoeff, long maxNegLag, long maxPosLag) {
	long ccLag = 0;
	corrCoeff = 0;

	for (long i = maxNegLag; i <= maxPosLag; i++) {
		double currCorrCoeff = SignalProcessing::getPearsonCorrCoeff(signalA, signalB, i);
		if (abs(currCorrCoeff) > abs(corrCoeff)) {
			corrCoeff = currCorrCoeff;
			ccLag = i;
		}
	}

	return ccLag;
}

long CorrelationFeatureSelection::getMaxCrossCorrelationRho(std::vector<double> &signalA, std::vector<double> &signalB, long startSample, long length, double& corrCoeff, long maxNegLag, long maxPosLag) {
	long ccLag = 0;
	double currCorrCoeff = 0;
	corrCoeff = 0;
	spearRslt rslt;

	//Look for the best correlation coefficient using a time lag from -500ms to 500ms with 0.5ms steps
	for (long i = maxNegLag; i <= maxPosLag; i++) {
		//currCorrCoeff = Correlation_c::getCorrCoeff(signalA, signalB, i);
		rslt = Correlation_c::getSpearmanCorrCoeff(signalA, signalB, i);
		currCorrCoeff = rslt.rs;
		double significance = rslt.probrs;

		if (abs(currCorrCoeff)> abs(corrCoeff)) {
			corrCoeff = currCorrCoeff;
			ccLag = i;
		}
	}

	return ccLag;
}

bool CorrelationFeatureSelection::relabelDataHFO(bool getCFS_Ripple, bool getCFS_FastRipple, bool getCFS_AllHFO, EOI_Features &epochFeatures) {
	unsigned nrEpochs = epochFeatures.m_time.size();
	unsigned nrEvents = 0;
	double negLabelVal = 0;

	/*std::string filename = m_patientPath + "/TrainingFiles/CFS/Labels";
	if (getCFS_Ripple) {
		filename += "_Ripples.txt";
	}
	else if (getCFS_FastRipple) {
		filename += "_FastRipples.txt";
	}
	remove(filename.c_str());
	std::ofstream labelsFile;
	labelsFile.open(filename, std::ios::app);    // open file for appending
	labelsFile << "Chann\t" << "Label\n";*/

	//Reset label for neighboing channels
	for (int channelIdx = 0; channelIdx < epochFeatures.m_EOI_Mark.size(); ++channelIdx) {
//#pragma omp parallel for reduction(+:nrEvents)
		for (int epochIdx = 0; epochIdx < nrEpochs; ++epochIdx) {

			double label = epochFeatures.m_EOI_Mark[channelIdx][epochIdx];
			//labelsFile << channelIdx << "\t" << label << "\n";

			bool is_R_FR = (compareDoubles(label, RIPPLE_FR_LABEL) || compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL));
			bool isRipple = (compareDoubles(label, RIPPLE_LABEL) || compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL));
			bool isFastRipple = (compareDoubles(label, FAST_RIPPLE_LABEL) || compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL));
			bool isHFO = is_R_FR || isRipple || isFastRipple;

			bool is_R_FR_Neighb = m_multichannel && (compareDoubles(label, RIPPLE_FR_LABEL + PROPAGATED_ANNOT_LABEL) || compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL));
			bool isRippleNeighb = m_multichannel && (compareDoubles(label, RIPPLE_LABEL + PROPAGATED_ANNOT_LABEL) || compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL));
			bool isFastRippleNeighb = m_multichannel && (compareDoubles(label, FAST_RIPPLE_LABEL + PROPAGATED_ANNOT_LABEL) || compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL));
			bool isHFO_Neighb = (is_R_FR_Neighb || isRippleNeighb || isFastRippleNeighb);

			if (getCFS_Ripple) {
				if (isRipple || isRippleNeighb) {
					label = 1;
					nrEvents++;
				}
				else {
					label = negLabelVal;
				}
			}
			else if (getCFS_FastRipple) {
				if (isFastRipple || isFastRippleNeighb) {
					label = 1;
					nrEvents++;
				}
				else {
					label = negLabelVal;
				}
			}
			else if (getCFS_AllHFO) {
				if (isHFO || isHFO_Neighb) {
					label = 1;
					nrEvents++;
				}
				else {
					label = negLabelVal;
				}
			}

			epochFeatures.m_EOI_Mark[channelIdx][epochIdx] = label;
		}
	}

	//labelsFile.close();
	return nrEvents > 0;
}

bool CorrelationFeatureSelection::relabelDataSpikesHFO(EOI_Features &epochFeatures) {
	unsigned nrEpochs = epochFeatures.m_time.size();
	unsigned nrEvents = 0;
	double negLabelVal = 0;

	//Reset label for neighboing channels
	for (int channelIdx = 0; channelIdx < m_nrChannels; ++channelIdx) {
#pragma omp parallel for reduction(+:nrEvents)
		for (int epochIdx = 0; epochIdx < nrEpochs; ++epochIdx) {

			double label = epochFeatures.m_EOI_Mark[channelIdx][epochIdx];
			bool is_Spike_R_FR = compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL);
			bool is_Spike_Ripple =  compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL);
			bool is_Spike_FastRipple = compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL);
			bool isSpikeHFO = (is_Spike_R_FR || is_Spike_Ripple || is_Spike_FastRipple);

			bool is_Spike_R_FR_Neighb = m_multichannel && compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool is_Spike_Ripple_Neighb = m_multichannel && compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool is_Spike_FastRipple_Neighb = m_multichannel && compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool isSpikeHFO_Neighb = (is_Spike_R_FR_Neighb || is_Spike_Ripple_Neighb || is_Spike_FastRipple_Neighb);

			if (isSpikeHFO || isSpikeHFO_Neighb) {
				label = 1;
				nrEvents++;
			}
			else {
				label = negLabelVal;
			}
			epochFeatures.m_EOI_Mark[channelIdx][epochIdx] = label;
		}
	}

	return nrEvents > 0;
}

bool CorrelationFeatureSelection::relabelDataSpikesAlone(EOI_Features &epochFeatures) {
	unsigned nrEpochs = epochFeatures.m_time.size();
	unsigned nrEvents = 0;
	double negLabelVal = 0;

	//Reset label for neighboing channels
	for (int channelIdx = 0; channelIdx < m_nrChannels; ++channelIdx) {
#pragma omp parallel for reduction(+:nrEvents)
		for (int epochIdx = 0; epochIdx < nrEpochs; ++epochIdx) {

			double label = epochFeatures.m_EOI_Mark[channelIdx][epochIdx];
			bool isSpikeAlone = compareDoubles(label, SPIKE_LABEL);
			bool isSpikeAlone_Neighb = m_multichannel && compareDoubles(label, SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);

			if (isSpikeAlone || isSpikeAlone_Neighb) {
				label = 1;
				nrEvents++;
			}
			else {
				label = negLabelVal;
			}
			epochFeatures.m_EOI_Mark[channelIdx][epochIdx] = label;
		}
	}

	return nrEvents > 0;
}

bool CorrelationFeatureSelection::relabelDataAllSpikes(EOI_Features &epochFeatures) {
	unsigned nrEpochs = epochFeatures.m_time.size();
	unsigned nrEvents = 0;
	double negLabelVal = 0;

	//Reset label for neighboing channels
	for (int channelIdx = 0; channelIdx < m_nrChannels; ++channelIdx) {
#pragma omp parallel for reduction(+:nrEvents)
		for (int epochIdx = 0; epochIdx < nrEpochs; ++epochIdx) {

			double label = epochFeatures.m_EOI_Mark[channelIdx][epochIdx];
			
			bool isSpikeAlone = compareDoubles(label, SPIKE_LABEL);
			bool isSpikeAlone_Neighb = m_multichannel && compareDoubles(label, SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);

			bool is_Spike_R_FR = compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL);
			bool is_Spike_Ripple = compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL);
			bool is_Spike_FastRipple = compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL);
			bool isSpikeHFO = (is_Spike_R_FR || is_Spike_Ripple || is_Spike_FastRipple);

			bool is_Spike_R_FR_Neighb = m_multichannel && compareDoubles(label, RIPPLE_FR_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool is_Spike_Ripple_Neighb = m_multichannel && compareDoubles(label, RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool is_Spike_FastRipple_Neighb = m_multichannel && compareDoubles(label, FAST_RIPPLE_LABEL + SPIKE_LABEL + PROPAGATED_ANNOT_LABEL);
			bool isSpikeHFO_Neighb = is_Spike_R_FR_Neighb || is_Spike_Ripple_Neighb || is_Spike_FastRipple_Neighb;

			bool isSpikeAllTypes = (isSpikeAlone || isSpikeHFO);
			bool isSpikeAllTypes_Neighb = (isSpikeAlone_Neighb || isSpikeHFO_Neighb);

			if (isSpikeAllTypes || isSpikeAllTypes_Neighb) {
				label = 1;
				nrEvents++;
			}
			else {
				label = negLabelVal;
			}
			epochFeatures.m_EOI_Mark[channelIdx][epochIdx] = label;
		}
	}

	return nrEvents > 0;
}

bool CorrelationFeatureSelection::vectorizeChannels(EOI_Features &epochFeatures, EOI_Features &vectorizedFeatureChannels) {

	unsigned nrFeatures = epochFeatures.features.size();
	unsigned nrEpochs = epochFeatures.m_time.size();

	vectorizedFeatureChannels.setMatrices(1);
	int channZero = 0;
	int totNrPosEvents = 0;

	for (long featIdx = 0; featIdx < nrFeatures; ++featIdx) {
		for (int channIdx = 0; channIdx < m_nrChannels; ++channIdx) {
			int nrPosEvents = 0;
			for (long sampleIdx = 0; sampleIdx < nrEpochs; ++sampleIdx) {
				double epochLabel = epochFeatures.m_EOI_Mark[channIdx][sampleIdx];
				if (epochLabel > 0)
					nrPosEvents++;
					totNrPosEvents++;
			}
			if (nrPosEvents >= 10) {
				vectorizedFeatureChannels.features.at(featIdx)->at(channZero).insert(std::end(vectorizedFeatureChannels.features.at(featIdx)->at(channZero)), std::begin(epochFeatures.features.at(featIdx)->at(channIdx)), std::end(epochFeatures.features.at(featIdx)->at(channIdx)));
				vectorizedFeatureChannels.m_EOI_Mark.at(channZero).insert(std::end(vectorizedFeatureChannels.m_EOI_Mark.at(channZero)), std::begin(epochFeatures.m_EOI_Mark.at(channIdx)), std::end(epochFeatures.m_EOI_Mark.at(channIdx)));
				if (vectorizedFeatureChannels.m_time.empty()) {
					vectorizedFeatureChannels.m_time.insert(std::end(vectorizedFeatureChannels.m_time), std::begin(epochFeatures.m_time), std::end(epochFeatures.m_time));
				}
			}
		}
	}

	return totNrPosEvents > 10;
}