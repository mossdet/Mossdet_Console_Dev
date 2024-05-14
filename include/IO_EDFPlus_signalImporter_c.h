// 26.09.2016 created by Daniel Lachner

#ifndef __IO_EDFPlus_SignalImporter_c_H__
#define __IO_EDFPlus_SignalImporter_c_H__

#include <vector>
#include <omp.h>
#include <sstream>

#include "DataTypes.h"

class  IO_EDFPlus_SignalImporter_c  
{
public:
	IO_EDFPlus_SignalImporter_c(const std::string strFilename, std::vector<std::string > selectedUnipolarChanns, std::vector<std::string > selectedBipolarChanns);
    virtual ~IO_EDFPlus_SignalImporter_c();

public:
	// public method to get data
	bool ReadSamples(long long lStartSample, long long lNumberSamples, matrixStd& signalSamples);
	bool ReadReformattedSamples(long long lStartSample, long long lNumberSamples, matrixStd& reformattedSamples);

	bool WriteReformattedSamples(long long lStartTime, long long lEndTime);

	//virtual bool setSelectedChannels(std::vector<unsigned> selection);

	// start date
	unsigned getChannelLength(unsigned ch);
	unsigned getChannelSamplingRate(unsigned ch);

	bool isValidImporter();

	double getSamplingRate() const;
	bool setSamplingRate(double samplingRate);

	unsigned long getSampleCount() const;
	bool setSampleCount(_int64 lSampleCount);

	unsigned int getNumberUnipolarChannels() const;
	unsigned int getNumberReformattedMontages() const;

	std::string getFileName();

	void setSelectedUnipolarChannels(std::vector<std::string > selectedChanns) {
		m_selectedUnipolarChannels = selectedChanns;
	}
	void setSelectedBipolarChannels(std::vector<std::string > selectedChanns) {
		m_selectedBipolarChannels = selectedChanns;
	}

	bool readUnipolarLabels(std::vector<ContactNames> &unipolarContactNames);
	bool readMontageLabels(std::vector<MontageNames> &montageLabels);

	DateAndTime getDateAndTime() {return m_fileStartDateAndTime;}
	std::string getPatientNameAndPath(void);

protected: 
	bool readHeader();
	inline _int64 getBlockLength() { return m_nBlockLength; };
	inline bool setBlockLength(_int64 nBlockLength) { m_nBlockLength = nBlockLength; return true; };

protected:
	_int64			m_nBlockLength;	// Block lenght in seconds
	long			m_lMinDig;			// minimum digital value
	long			m_lMaxDig;			// maximum digital value
	double			m_dMinPhys;			// minimum physical value
	double			m_dMaxPhys;			// maximum physical value
	double			m_dSamplingRate;	// signal sampling rate
	_int64	m_lSampleCount;		// sample count
	std::string m_strFilename;			// file name of signal file 
	std::string m_patientName;
	std::string m_patientPath;

	bool m_bValidImporter;
	std::vector<ContactNames> m_unipolarContacts;
	std::vector<MontageNames> m_montages;
	std::vector<std::string> m_selectedUnipolarChannels;
	std::vector<std::string> m_selectedBipolarChannels;
	std::vector<std::string> m_readChannelLabels;

	DateAndTime m_fileStartDateAndTime;

};

#endif // __IOAbstractSignalImporter_c__H__
