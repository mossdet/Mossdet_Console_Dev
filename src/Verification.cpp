#pragma once

#include <string>       // std::string
#include <sstream>      // std::stringstream
#include <fstream>
#include <random>
#include <map>
#include <direct.h>
#include <ctime>
#include <cstdlib>

#include <windows.h>

#include "Verification.h"

using namespace System::Runtime::InteropServices;

// random generator function:
int myrandom(int i) { return std::rand() % i; }

int compareEvent(EOI_Event a, EOI_Event b) {
	return a.startTime < b.startTime;
}

Verification::Verification(System::String^ strInputFileName, System::String^ strDirName, int nrSamples, double ci, double me) {

	m_strFileName = (char*)(void*)Marshal::StringToHGlobalAnsi(strInputFileName);
	m_outputDirectory = (char*)(void*)Marshal::StringToHGlobalAnsi(strDirName);
	getPatientNameAndPath();

	m_nrSamples = nrSamples;

	m_confInterval = ci;
	m_marginOfError = me;

	m_dEarliestDetectTime = 0;
	m_dLatestDetectTime = 0;

}
Verification::~Verification() {

}

bool Verification::getRandomSamples() {
	
	double dSignalLength = 0;
	double eventTypeDuration = 0;

	std::srand(unsigned(std::time(0)));

	bool frDetections, rippleDetections, hfoDetections, spikeDetections;

	if (m_strFileName.find("ast") != std::string::npos) {
		eventTypeDuration = 0.025;
		frDetections = true;
	}
	else if ((m_strFileName.find("ipple") != std::string::npos) && !frDetections) {
		eventTypeDuration = 0.05;
		rippleDetections = true;
	}
	else if (m_strFileName.find("HFO") != std::string::npos) {
		eventTypeDuration = 0.05;
		hfoDetections = true;
	}
	else if (m_strFileName.find("ike") != std::string::npos) {
		eventTypeDuration = 0.100;
		spikeDetections = true;
	}

	std::ifstream detections;
	detections.open(m_strFileName.c_str());
	long nrPosEvents = 0;
	//read channels, file start and end times
	if (detections) {
		detections.precision(16);
		std::string s_length, s_header;
		std::getline(detections, s_length);	// signal duration
		std::getline(detections, s_header);	// header
		while (!detections.eof()) {
			std::string line;
			std::string dataLabel;
			std::getline(detections, line);
			std::stringstream ss(line);

			ss >> dataLabel;
			if (dataLabel.compare("d:") == 0) {

				EOI_Event eoiMark;
				std::string startTimeDay, endTimeDay;
				int dayNr, match;

				ss >> eoiMark.description;
				ss >> eoiMark.channelName;
				ss >> eoiMark.channelNr;
				ss >> eoiMark.startTime;
				ss >> eoiMark.endTime;
				ss >> eoiMark.duration;
				ss >> dayNr;
				ss >> startTimeDay;
				ss >> endTimeDay;
				ss >> match;

				if (!eoiMark.channelName.empty() && eoiMark.startTime > 0.0 && eoiMark.endTime > 0.0) {
					nrPosEvents++;
					m_channels[eoiMark.channelNr] = eoiMark.channelName;
					if (nrPosEvents == 1) {
						m_dEarliestDetectTime = eoiMark.startTime;
					}
					else {
						if (eoiMark.startTime < m_dEarliestDetectTime) {
							m_dEarliestDetectTime = eoiMark.startTime;
						}
						if (eoiMark.endTime > m_dLatestDetectTime) {
							m_dLatestDetectTime = eoiMark.endTime;
						}

					}

				}
			}
		}
		detections.close();


		// Save the positive events
		m_allPositiveEvents.resize(m_channels.size());
		m_allNegativeEvents.resize(m_channels.size());
		detections.open(m_strFileName.c_str());
		detections.precision(16);
		detections.seekg(0, detections.beg);		
		std::getline(detections, s_length);	// signal duration
		std::getline(detections, s_header);	// header
		while (!detections.eof()) {
			std::string line;
			std::string dataLabel;
			std::getline(detections, line);
			std::stringstream ss(line);

			ss >> dataLabel;
			if (dataLabel.compare("d:") == 0) {

				EOI_Event eoiMark;
				std::string startTimeDay, endTimeDay;
				int dayNr, match;

				ss >> eoiMark.description;
				ss >> eoiMark.channelName;
				ss >> eoiMark.channelNr;
				ss >> eoiMark.startTime;
				ss >> eoiMark.endTime;
				ss >> eoiMark.duration;
				ss >> dayNr;
				ss >> startTimeDay;
				ss >> endTimeDay;
				ss >> match;

				if (!eoiMark.channelName.empty() && eoiMark.startTime > 0.0 && eoiMark.endTime > 0.0) {
					m_allPositiveEvents[eoiMark.channelNr].push_back(eoiMark);
				}
			}
		}
		detections.close();
	}
	

	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> randChannGen(0, m_channels.size() - 1); // define the range
	std::uniform_int_distribution<> randTimeGen(m_dEarliestDetectTime, m_dLatestDetectTime - 1); // define the range
	long nrNegativeEvents = 0;
	while (nrNegativeEvents < nrPosEvents) {

		double randomStartTime = randTimeGen(eng) + (double)(std::rand() % 1000) / 1000.0;
		int randomChannel = randChannGen(eng);

		EOI_Event eoiNegativeMark;
		eoiNegativeMark.description = "NegEvent";
		eoiNegativeMark.startTime = randomStartTime;
		eoiNegativeMark.endTime = eoiNegativeMark.startTime + eventTypeDuration;
		eoiNegativeMark.duration = eventTypeDuration;
		eoiNegativeMark.channelNr = randomChannel;
		eoiNegativeMark.channelName = m_channels.at(randomChannel);

		for (long posEvIdx = 0; posEvIdx < m_allPositiveEvents[randomChannel].size(); ++posEvIdx) {
			
			double posEvStart = m_allPositiveEvents[randomChannel][posEvIdx].startTime;
			double posEvEnd = m_allPositiveEvents[randomChannel][posEvIdx].endTime;
			double posEvDuration = m_allPositiveEvents[randomChannel][posEvIdx].duration;

			double shared = 0;
			if (posEvStart <= eoiNegativeMark.startTime && posEvEnd >= eoiNegativeMark.startTime) {		//A
				if (posEvEnd <= eoiNegativeMark.endTime)
					shared = posEvEnd - eoiNegativeMark.startTime;
				else
					shared = eoiNegativeMark.endTime - eoiNegativeMark.startTime;
			}
			else if (posEvStart >= eoiNegativeMark.startTime && posEvStart <= eoiNegativeMark.endTime) {	//B
				if (posEvEnd >= eoiNegativeMark.endTime)
					shared = eoiNegativeMark.endTime - posEvStart;
				else
					shared = posEvEnd - posEvStart;
			}

			if (shared < 0.0000001) {
				m_allNegativeEvents[randomChannel].push_back(eoiNegativeMark);
				nrNegativeEvents++;
				break;
			}
		}
	}

	for (int chIdx = 0;chIdx < m_channels.size(); ++chIdx) {
		std::random_shuffle(m_allNegativeEvents[chIdx].begin(), m_allNegativeEvents[chIdx].end(), myrandom);
		std::random_shuffle(m_allPositiveEvents[chIdx].begin(), m_allPositiveEvents[chIdx].end(), myrandom);
	}

	return true;
}

bool Verification::generateVerificationFile(void) {

	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> randChannGen(0, m_channels.size() - 1); // define the range

	std::vector<EOI_Event> verificationSamples;

	long classSampIdx = 0;
	for (long sampIdx = 0; sampIdx < m_nrSamples/2; ++sampIdx) {
		int randChann = randChannGen(eng);
		std::uniform_int_distribution<> randPosSampGen(0, m_allPositiveEvents[randChann].size() - 1); // define the range
		long randPosSamp = randPosSampGen(eng);
		std::uniform_int_distribution<> randNegSampGen(0, m_allNegativeEvents[randChann].size() - 1); // define the range
		long randNegSamp = randNegSampGen(eng);
		verificationSamples.push_back(m_allPositiveEvents[randChann][randPosSamp]);
		verificationSamples.push_back(m_allNegativeEvents[randChann][randNegSamp]);
	}

	std::sort(verificationSamples.begin(), verificationSamples.end(), compareEvent);


	std::string outputPath = m_outputDirectory + "/MOSSDET_Output";
	mkdir(outputPath.c_str());
	m_confInterval;
	m_marginOfError;

	std::string filename = outputPath + "/" + m_patientName + "_CI_" + std::to_string((int)m_confInterval) + "%_ME_" + std::to_string((int)m_marginOfError) + "%_detectionsToVerify.txt";
	std::ofstream verifyDetectionsFile;
	remove(filename.c_str());
	verifyDetectionsFile.open(filename);
	verifyDetectionsFile << "Data\t" << "Description\t" << "ChannelName\t" << "ChannelNr.\t" << "StartTime(s)\t" << "EndTime(s)\t" << "Duration(s)\t" << "DayNr.\t" << "StartTime(day)\t" << "EndTime(day)\t" << "Match\n";

	for (long di = 0; di < verificationSamples.size(); ++di) {
		verifyDetectionsFile << "d:" << "\t"
			<< verificationSamples[di].description << "\t"
			<< verificationSamples[di].channelName << "\t"
			<< verificationSamples[di].channelNr << "\t"
			<< verificationSamples[di].startTime << "\t"
			<< verificationSamples[di].endTime << "\t"
			<< verificationSamples[di].duration << "\t"
			<< 0 << "\t"
			<< 0 << "\t"
			<< 0 << "\t"
			<< 0 << "\n";
	}
	
	verifyDetectionsFile.close();
	return true;
}


std::string Verification::getPatientNameAndPath() {
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