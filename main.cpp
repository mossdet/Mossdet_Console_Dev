#include <string>
#include <vector>
#include <direct.h>
#include <chrono>
#include <fstream>
#include <iostream>

#include "HFO_Detector.h"
#include "SleepSpindleDetector.h"
#include "SignalProcessing.h"
#include "IO_EDFPlus_SignalImporter_c.h"
#include "MatlabFilesRW.h"

#include <omp.h>

//#define NO_INPUT_ARGS

std::string m_strOutputDirectory;

int main(int argc, char **argv) {

#ifndef	NO_INPUT_ARGS
	if (argc == 10) {
#endif
		// filename, output path, 
		auto start_time = std::chrono::high_resolution_clock::now();

		std::string strInputFileName;
		std::string detectorPath;
		double startTime;
		double endTime;
		double samplingRate;
		std::string eoiName;
		bool verbose = false;
		bool saveTextFile = false;

#ifndef	NO_INPUT_ARGS
		strInputFileName = argv[1];
		detectorPath = argv[2];
		m_strOutputDirectory = argv[3];
		startTime = atoi(argv[4]);
		endTime = atoi(argv[5]);
		samplingRate = atoi(argv[6]);
		eoiName = argv[7];
		verbose = atoi(argv[8]) > 0;
		saveTextFile = atoi(argv[9]) > 0;

#else
		strInputFileName = "F:\\MOSSDET_Development\\MOSSDET_c\\mattest.mat";
		detectorPath = "F:\\MOSSDET_Development\\MOSSDET_c";
		m_strOutputDirectory = "F:\\MOSSDET_Development\\MOSSDET_c";
		startTime = 0;
		endTime = 3600; 60 * 60 * 24 * 365;
		samplingRate = 2000;
		eoiName = "SleepSpindles"; //"HFO+IES", "SleepSpindles"
		//eoiName = "HFO+IES"; //"HFO+IES", "SleepSpindles"
		verbose = true;
		saveTextFile = true;

#endif

		if (verbose) {
			std::cout << "MOSSDET Detector" << std::endl;
			std::cout << "InputParameters: DataFile\t SVM_Path\t OutputPath\t StartTime\t EndTime\t SamplingRate\t EOI_Type\t" << std::endl << std::endl;
			std::cout << "Possible EOI_Types: HFO+IES\t SleepSpindles\n";
			std::cout << "HFO+IES : Detects all Ripples, FastRipples and Spikes\n";
			std::cout << "SleepSpindles : Detects only sleep-spindles\n\n";

#ifndef	NO_INPUT_ARGS
			std::cout << "argc:" << argc << std::endl;
			std::cout << "Filename:" << argv[1] << std::endl;
			std::cout << "DetectorPath:" << argv[2] << std::endl;
			std::cout << "OutPath:" << argv[3] << std::endl;
			std::cout << "StarTime:" << atoi(argv[4]) << std::endl;
			std::cout << "EndTime:" << atoi(argv[5]) << std::endl;
			std::cout << "SamplingRate:" << atoi(argv[6]) << std::endl;
			std::cout << "EOI:" << argv[7] << std::endl << std::endl;
			std::cout << "Verbose:" << argv[8] << std::endl << std::endl;
			std::cout << "SaveTextFile:" << argv[9] << std::endl << std::endl;
#endif
		}

		int output = -1;
		if (eoiName.compare("HFO+IES")== 0) {
			std::cout << "Detect HFO+IES\n";
			std::vector<std::string> selectedUnipolarChanels, selectedBipolarChanels;
			HFO_Detector hfoRippleDetector(strInputFileName, detectorPath, m_strOutputDirectory, 1, false, &selectedUnipolarChanels, &selectedBipolarChanels, eoiName, false, startTime, endTime, false, samplingRate, verbose, saveTextFile);
			output = hfoRippleDetector.characterizeEEG();
		}
		else if (eoiName.compare("SleepSpindles") == 0) {
			std::cout << "SleepSpindles\n";
			std::vector<std::string> selectedUnipolarChanels, selectedBipolarChanels;
			SleepSpindleDetector ssDetector(strInputFileName, detectorPath, m_strOutputDirectory, 1, false, &selectedUnipolarChanels, &selectedBipolarChanels, eoiName, false, startTime, endTime, false, samplingRate, verbose, saveTextFile);
			//output = ssDetector.readEEG();
			output = ssDetector.characterizeEEG();
		}
		else {
			std::cout << "Wrong EOI Name!" << std::endl << std::endl;
			return 0;
		}


		if (output == 1) {
			std::cout << "OK" << std::endl;
		}
		else if (output == -1) {
			std::cout << "E#1" << std::endl;
		}
		else if (output == -2) {
			std::cout << "E#2" << std::endl;
		}

#ifndef	NO_INPUT_ARGS
	}
	else {
		printf("Wrong Input Parameters\n");
	}
#endif


	return 1;
}
