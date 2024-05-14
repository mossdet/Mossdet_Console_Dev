#pragma once

#include <string>
#include <math.h>
#include <map>

#include "DataTypes.h"

class Verification
{
public:
	Verification(System::String^ strInputFileName, System::String^ strDirName, int m_nrSamples, double ci, double me);
	~Verification();

	bool generateVerificationFile(void);
	std::string getPatientNameAndPath();
	bool getRandomSamples();

private:


	std::string m_strFileName;
	std::string m_outputDirectory;
	std::string m_patientName;
	std::string m_patientPath;
	int m_nrSamples;
	std::map<int, std::string> m_channels;
	double m_dEarliestDetectTime;
	double m_dLatestDetectTime;

	double m_confInterval;
	double m_marginOfError;

	std::vector<std::vector<EOI_Event>> m_allPositiveEvents;
	std::vector<std::vector<EOI_Event>> m_allNegativeEvents;
};