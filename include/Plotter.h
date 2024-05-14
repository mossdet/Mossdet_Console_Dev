#pragma once

#include <string>
#include <math.h>

#include "DataTypes.h"

class Plotter {
public:
	Plotter() {};
	~Plotter() {};

	bool Plot_3D_Profile(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::string &wdwTitle, double xtics, double xmin, double xmax, double ytics, double ymin, double ymax);
	static void plotXY(std::vector<double> &x, std::vector<double> &y, std::string &wdwTitle, double xtics, double ytics);
	static void writeVectorsToFile(std::vector<double> &x, std::vector<double> &y, std::string &filename);


private:
	std::string m_savedFunctionName;
	std::string m_eoiName;

	// Characterized Signal parameters
	std::string m_patientName;
	std::string m_patientPath;
	std::string m_filename;
	
	std::vector<EOI_Epoch> m_inputEvents;
	std::vector<std::vector<bool>> m_chSpecifificDetecEpochs;
	std::vector<EOI_Event> m_detectedEOI;
	std::vector<EOI_Event> m_visualEOI;
	std::vector<EOI_Event> m_visualEOI_Joined;

	//Initialize File Paths
	std::string m_outputPath;
	std::string m_patPath;
	std::string m_patDetectionFilesPath;
};