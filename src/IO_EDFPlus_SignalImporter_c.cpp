// 26.09.2016 created by Daniel Lachner-Piza

#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;

#include "IO_EDFPlus_SignalImporter_c.h"
#include "edflib.h"
# include <omp.h>


//////////////////////////////////////////////////////////////////////
// constructor, destructor

IO_EDFPlus_SignalImporter_c::IO_EDFPlus_SignalImporter_c(const std::string strFilename, std::vector<std::string > selectedUnipolarChanns, std::vector<std::string > selectedBipolarChanns)
{
	// read header and fill necesarry members
	m_strFilename = strFilename;
	m_selectedUnipolarChannels = selectedUnipolarChanns;
	m_selectedBipolarChannels = selectedBipolarChanns;

	m_bValidImporter = readHeader();

	getPatientNameAndPath();
}


IO_EDFPlus_SignalImporter_c::~IO_EDFPlus_SignalImporter_c()
{
}


// read the date
bool IO_EDFPlus_SignalImporter_c::ReadSamples(long long lStartSample, long long lNumberSamples, matrixStd& signalSamples)
{
	// check for valid importer
	if (!isValidImporter())
		return false;

	// check requested samples against available samples
	if (lStartSample + lNumberSamples > getSampleCount())
		lNumberSamples = getSampleCount() - lStartSample;

	std::string edfFile = getFileName();
	edf_hdr_struct edfhdr;
	_int64 headerReadOK = edfopen_file_readonly(edfFile.c_str(), &edfhdr, EDFLIB_DO_NOT_READ_ANNOTATIONS);
	edfhdr.filetype = EDFLIB_FILETYPE_EDF;

	if (headerReadOK == -1)
		return false;

	if (signalSamples.empty()) {				//First read, load all selected channels to memory
		
		unsigned int nSignals = !m_unipolarContacts.empty() ? m_unipolarContacts.size() : edfhdr.edfsignals;
		signalSamples.resize(nSignals);

		!m_unipolarContacts.empty() ? signalSamples.resize(m_unipolarContacts.size()) : signalSamples.resize(nSignals);
//#pragma omp parallel for
		for (int signalNr = 0; signalNr < nSignals; ++signalNr) {

			std::string chName = edfhdr.signalparam[m_unipolarContacts[signalNr].contactGlobalIdx].label;
			unsigned samplingFreq = edfhdr.signalparam[m_unipolarContacts[signalNr].contactGlobalIdx].smp_in_datarecord;
			_int64 totalSamples = lNumberSamples; edfhdr.signalparam[m_unipolarContacts[signalNr].contactGlobalIdx].smp_in_file;
			std::string dimensions = edfhdr.signalparam[m_unipolarContacts[signalNr].contactGlobalIdx].physdimension;
			std::string preFilter = edfhdr.signalparam[m_unipolarContacts[signalNr].contactGlobalIdx].prefilter;
			double *buf = new double[totalSamples];

			edfseek((int)edfhdr.handle, (int)m_unipolarContacts[signalNr].contactGlobalIdx, lStartSample, EDFSEEK_SET);//set start Sample
			_int64 out = edfread_physical_samples((int)edfhdr.handle, (int)m_unipolarContacts[signalNr].contactGlobalIdx, (int)totalSamples, buf);

			double dScale = -1;// 1E-6;
			signalSamples[signalNr].resize(totalSamples);
#pragma omp parallel for
			for (long i = 0; i < totalSamples; ++i) {
				double val = buf[i] * dScale;
				signalSamples[signalNr][i] = val;
			}

			delete[] buf;
		}
	}

	edfclose_file(edfhdr.handle);

	/*bool writeChannels = false;
	if (writeChannels) {
	std::ofstream channelFile;
	std::string filename = chName + ".txt";
	remove(filename.c_str());
	channelFile.open(filename);
	channelFile.precision(64);
	for (unsigned i = 0; i < nrSamples; ++i) {
	double val = buf[i] * 1E-6;
	channelFile << val << "\n";
	}
	channelFile.close();
	}*/

	return true;
}

bool IO_EDFPlus_SignalImporter_c::ReadReformattedSamples(long long lStartSample, long long lNumberSamples, matrixStd& reformattedSamples)
{
	// check for valid importer
	if (!isValidImporter())
		return false;

	// check requested samples against available samples
	if (lStartSample + lNumberSamples > getSampleCount())
		lNumberSamples = getSampleCount();

	if (reformattedSamples.empty()) {				//First read, load all selected channels to memory

		int nrMontages = (int)m_montages.size();
		reformattedSamples.resize(nrMontages);
		long nrSamples = lNumberSamples;// getSampleCount();
		for (long c = 0; c < nrMontages; ++c)
			reformattedSamples[c].resize(nrSamples);

		for (int montageNr = 0; montageNr < nrMontages; ++montageNr) {

			std::string edfFile = getFileName();
			edf_hdr_struct edfhdr;
			_int64 headerReadOK = edfopen_file_readonly(edfFile.c_str(), &edfhdr, EDFLIB_DO_NOT_READ_ANNOTATIONS);
			if (headerReadOK == -1)
				return false;

			double *bufA = new double[lNumberSamples];
			double *bufB = new double[lNumberSamples];

			edfseek((int)edfhdr.handle, (int)m_montages[montageNr].firstContactGlobalIdx, lStartSample, EDFSEEK_SET);//set start Sample
			edfseek((int)edfhdr.handle, (int)m_montages[montageNr].secondContactGlobalIdx, lStartSample, EDFSEEK_SET);//set start Sample

			_int64 outA = edfread_physical_samples((int)edfhdr.handle, (int)m_montages[montageNr].firstContactGlobalIdx, (int)lNumberSamples, bufA);
			_int64 outB = edfread_physical_samples((int)edfhdr.handle, (int)m_montages[montageNr].secondContactGlobalIdx, (int)lNumberSamples, bufB);

			if (outA == -1 || outB == -1)
				return false;
			
			double dScale = -1; 1E-6;
			
//#pragma omp parallel for
			for (long i = 0; i < lNumberSamples; ++i) {
				double valA = bufA[i];
				double valB = bufB[i];
				double val = (double)(valA - valB);
				val *= dScale;
				reformattedSamples[montageNr][i] = val;
			}

			delete[] bufA;
			delete[] bufB;
			
			edfclose_file(edfhdr.handle);
		}
	}

	return true;
}

// read header
bool IO_EDFPlus_SignalImporter_c::readHeader()
{

	std::string edfFile = getFileName();
	edf_hdr_struct edfhdr;
	int headerReadOK = edfopen_file_readonly(edfFile.c_str(), &edfhdr, EDFLIB_DO_NOT_READ_ANNOTATIONS); //EDFLIB_READ_ANNOTATIONS, EDFLIB_READ_ALL_ANNOTATIONS; EDFLIB_DO_NOT_READ_ANNOTATIONS
	//edfhdr.filetype = EDFLIB_FILETYPE_EDFPLUS; EDFLIB_FILETYPE_EDF;

	unsigned signalNr = 6;
	std::string chName = edfhdr.signalparam[signalNr].label;
	unsigned nrSamples = edfhdr.signalparam[signalNr].smp_in_file;
	std::string dimensions = edfhdr.signalparam[signalNr].physdimension;
	std::string preFilter = edfhdr.signalparam[signalNr].prefilter;
	
	m_fileStartDateAndTime.day = edfhdr.startdate_day;
	m_fileStartDateAndTime.month = edfhdr.startdate_month;
	m_fileStartDateAndTime.year = edfhdr.startdate_year;

	m_fileStartDateAndTime.hours = edfhdr.starttime_hour;
	m_fileStartDateAndTime.minutes = edfhdr.starttime_minute;
	m_fileStartDateAndTime.seconds = edfhdr.starttime_second;
	m_fileStartDateAndTime.milliSeconds = edfhdr.starttime_subsecond;

	// get number of signals
	unsigned int nSignals = edfhdr.edfsignals;
	// get number of blocks
	_int64 nBlocks = edfhdr.datarecords_in_file;
	// get block duration in seconds
	_int64 nDuration = edfhdr.datarecord_duration / (10 * 1000 * 1000); // only int values possible ?
	setBlockLength(nDuration);

	// Get Sampling Rate and Count
	// compute sample count
	_int64 sampleCount = edfhdr.signalparam[1].smp_in_file;
	if (!setSampleCount(sampleCount))
		return false;
	// set sampling rate 
	double dataRecordDurationInSec = (edfhdr.datarecord_duration * 100.0) / (1000.0 * 1000.0 * 1000.0);
	double samplingRate = (double)edfhdr.signalparam[0].smp_in_datarecord / dataRecordDurationInSec;
	if (!setSamplingRate(samplingRate))
		return false;

	std::map<double, int> samplingRatesOcurrence;
	for (unsigned int n = 0; n < nSignals; n++) {
		double sr = edfhdr.signalparam[n].smp_in_datarecord;
		if (samplingRatesOcurrence.find(sr) == samplingRatesOcurrence.end())
			samplingRatesOcurrence[sr] = 1;
		else
			samplingRatesOcurrence[sr] += 1;
	}
	//Get Unipolar Channel Names
	for (unsigned int n = 0; n < nSignals; n++) {

		//clean read labels
		std::string strLabel = edfhdr.signalparam[n].label;
		std::size_t found = strLabel.find("EEG ");
		if (found != std::string::npos) {
			strLabel.erase(found, found + 4);
		}
		found = strLabel.find("EEG");
		if (found != std::string::npos) {
				strLabel.erase(found, found + 3);
		}
		strLabel.erase(std::remove_if(strLabel.begin(), strLabel.end(), ::isspace), strLabel.end());

		std::string electrodeName, contactNrString;

		std::vector<std::string> numbers;
		for (int nrIdx = 0; nrIdx < 10; nrIdx++) {
			numbers.push_back(std::to_string(nrIdx));
		}

		std::string::size_type minIdx = 1000000, maxIdx = 0, length = 0;
		for (int nrIdx = 0; nrIdx < 10; nrIdx++) {
			std::string::size_type fi = strLabel.find(numbers[nrIdx]);
			std::string::size_type si = strLabel.rfind(numbers[nrIdx]);

			if (fi != std::string::npos) {
				if (fi < minIdx) {
					minIdx = fi;
				}
			}

			if (si != std::string::npos) {
				if (si > maxIdx) {
					maxIdx = si;
				}
			}
		}

		int contactNr;
		if (minIdx != 1000000) {
			length = (maxIdx - minIdx) + 1;
			electrodeName = strLabel.substr(0, minIdx);
			contactNrString = strLabel.substr(minIdx, length);
			std::stringstream ss(contactNrString);
			ss >> contactNr;
		}
		else {
			electrodeName = strLabel;
			contactNr = 0;
		}

		ContactNames contact;
		contact.contactGlobalIdx = n;
		contact.contactNr = contactNr;
		contact.electrodeName = electrodeName;
		contact.contactName = strLabel;

		bool channelNotSelected = true;
		// if there are any selected channels, skip channels which are not selected
		if (!m_selectedUnipolarChannels.empty()) {
			std::string edfChannName = contact.contactName;
			for (int i = 0; i < m_selectedUnipolarChannels.size(); ++i) {
				std::string selectedChannName = m_selectedUnipolarChannels[i];
				selectedChannName.erase(std::remove_if(selectedChannName.begin(), selectedChannName.end(), ::isspace), selectedChannName.end());
				edfChannName.erase(std::remove_if(edfChannName.begin(), edfChannName.end(), ::isspace), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '-'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '-'), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '_'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '_'), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '#'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '#'), edfChannName.end());

				int compareResult = selectedChannName.compare(edfChannName);
				if (compareResult == 0) {
					channelNotSelected = false;
					break;
				}
			}
		}
		else {
			channelNotSelected = false;
		}
		if (channelNotSelected) {
			continue;
		}
		m_unipolarContacts.push_back(contact);
		/*}
		else
			bool stop = true;*/
	}

	edfclose_file(edfhdr.handle);

	int montageNr = 0;
	//read Montages' Names
	for (int channIdx = 0; channIdx < m_unipolarContacts.size()-1; ++channIdx) {
		MontageNames montage;
		montage.firstElectrodeName = m_unipolarContacts[channIdx].electrodeName;
		montage.firstContactNr = m_unipolarContacts[channIdx].contactNr;
		montage.firstContactGlobalIdx = m_unipolarContacts[channIdx].contactGlobalIdx;

		montage.secondElectrodeName = m_unipolarContacts[channIdx + 1].electrodeName;
		montage.secondContactNr = m_unipolarContacts[channIdx + 1].contactNr;
		montage.secondContactGlobalIdx = m_unipolarContacts[channIdx + 1].contactGlobalIdx;

		montage.montageName = m_unipolarContacts[channIdx].electrodeName + std::to_string(m_unipolarContacts[channIdx].contactNr) + "-" + m_unipolarContacts[channIdx+1].electrodeName + std::to_string(m_unipolarContacts[channIdx+1].contactNr);

		montage.montageMOSSDET_Nr = montageNr;

		// if there are any selected channels, skip channels which are not selected
		/*if (!m_selectedBipolarChannels.empty()) {
			if (std::find(m_selectedBipolarChannels.begin(), m_selectedBipolarChannels.end(), montage.montageName) == m_selectedBipolarChannels.end()) {
				continue;
			}
		}*/

		bool channelNotSelected = true;
		if (!m_selectedBipolarChannels.empty()) {
			std::string edfChannName = montage.montageName;
			for (int i = 0; i < m_selectedBipolarChannels.size(); ++i) {
				std::string selectedChannName = m_selectedBipolarChannels[i];
				selectedChannName.erase(std::remove_if(selectedChannName.begin(), selectedChannName.end(), ::isspace), selectedChannName.end());
				edfChannName.erase(std::remove_if(edfChannName.begin(), edfChannName.end(), ::isspace), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '-'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '-'), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '_'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '_'), edfChannName.end());

				selectedChannName.erase(std::remove(selectedChannName.begin(), selectedChannName.end(), '#'), selectedChannName.end());
				edfChannName.erase(std::remove(edfChannName.begin(), edfChannName.end(), '#'), edfChannName.end());

				int compareResult = selectedChannName.compare(edfChannName);
				if (compareResult == 0) {
					channelNotSelected = false;
					break;
				}
			}
		}
		else {
			channelNotSelected = false;
		}
		if (channelNotSelected) {
			continue;
		}
		if (montage.firstElectrodeName.compare(montage.secondElectrodeName) == 0) {
			m_montages.push_back(montage);
			montageNr++;
		}
		
	}

	// header is valid
	return headerReadOK == 0;
}

unsigned IO_EDFPlus_SignalImporter_c::getChannelLength(unsigned ch) {
	std::string edfFile = getFileName();
	edf_hdr_struct edfhdr;
	int headerReadOK = edfopen_file_readonly(edfFile.c_str(), &edfhdr, EDFLIB_READ_ALL_ANNOTATIONS);
	_int64 nrSamples = edfhdr.signalparam[ch].smp_in_file;
	edfclose_file(edfhdr.handle);

	return true;
}

unsigned IO_EDFPlus_SignalImporter_c::getChannelSamplingRate(unsigned ch) {
	std::string edfFile = getFileName();
	edf_hdr_struct edfhdr;
	int headerReadOK = edfopen_file_readonly(edfFile.c_str(), &edfhdr, EDFLIB_READ_ALL_ANNOTATIONS);
	unsigned samplingFreq = edfhdr.signalparam[ch].smp_in_datarecord/2;
	edfclose_file(edfhdr.handle);

	return samplingFreq;
}

// valid importer
bool IO_EDFPlus_SignalImporter_c::isValidImporter() {
	double samplingRate = getSamplingRate();
	unsigned sampleCount = getSampleCount();
	unsigned nrChs = getNumberUnipolarChannels();
	// check sampling rate, sample count and number channels
	if (getSamplingRate() <= 0.0 || getSampleCount() == 0 || getNumberUnipolarChannels() == 0)
		m_bValidImporter = false;

	return m_bValidImporter;
};

double IO_EDFPlus_SignalImporter_c::getSamplingRate() const{
		return m_dSamplingRate;
}

bool IO_EDFPlus_SignalImporter_c::setSamplingRate(double dSamplingRate) {
	m_dSamplingRate = dSamplingRate;

	if (m_dSamplingRate < 0) {
		m_bValidImporter = false;
		return false;
	}

	return true;
}

// sample count
unsigned long IO_EDFPlus_SignalImporter_c::getSampleCount() const {
	return (unsigned long)m_lSampleCount;
}

bool IO_EDFPlus_SignalImporter_c::setSampleCount(_int64 lSampleCount) {
	m_lSampleCount = lSampleCount;

	if (m_lSampleCount == 0) {
		m_bValidImporter = false;
		return false;
	}

	return true;
}

// number of channels

unsigned int IO_EDFPlus_SignalImporter_c::getNumberUnipolarChannels() const {
	return (unsigned int)m_unipolarContacts.size();
}

unsigned int IO_EDFPlus_SignalImporter_c::getNumberReformattedMontages() const {
	return (unsigned int)m_montages.size();
}

// file name
std::string IO_EDFPlus_SignalImporter_c::getFileName() {
	return m_strFilename;
}



bool IO_EDFPlus_SignalImporter_c::readUnipolarLabels(std::vector<ContactNames> &unipolarContactNames) {
	if (m_unipolarContacts.size() == 0)
		return false;

	for (int i = 0; i < m_unipolarContacts.size(); ++i)
		unipolarContactNames.push_back(m_unipolarContacts[i]);

	return true;
}

bool IO_EDFPlus_SignalImporter_c::readMontageLabels(std::vector<MontageNames> &montageLabels) {
	if (m_montages.size() == 0)
		return false;

	for (int i = 0; i < m_montages.size(); ++i)
		montageLabels.push_back(m_montages[i]);

	return true;
}


bool IO_EDFPlus_SignalImporter_c::WriteReformattedSamples(long long lStartTime, long long lEndTime)
{
	long long lStartSampleRead = lStartTime * getSamplingRate();
	long long lNumberSamplesRead = lEndTime * getSamplingRate();

	matrixStd reformattedSamples;

	// check for valid importer
	if (!isValidImporter())
		return false;

	// check requested samples against available samples
	if (lStartSampleRead + lNumberSamplesRead > getSampleCount())
		lNumberSamplesRead = getSampleCount();

	if (reformattedSamples.empty()) {				//First read, load all selected channels to memory
		int nrMontages = (int)m_montages.size();
		reformattedSamples.resize(nrMontages);
		long nrSamples = lNumberSamplesRead;// getSampleCount();
		for (long c = 0; c < nrMontages; ++c)
			reformattedSamples[c].resize(nrSamples);

		for (int montageNr = 0; montageNr < nrMontages; ++montageNr) {
			edf_hdr_struct edfhdr;
			_int64 headerReadOK = edfopen_file_readonly(getFileName().c_str(), &edfhdr, EDFLIB_DO_NOT_READ_ANNOTATIONS);
			if (headerReadOK == -1)
				return false;

			double *bufA = new double[lNumberSamplesRead];
			double *bufB = new double[lNumberSamplesRead];
			edfseek((int)edfhdr.handle, (int)m_montages[montageNr].firstContactGlobalIdx, lStartSampleRead, EDFSEEK_SET);//set start Sample
			edfseek((int)edfhdr.handle, (int)m_montages[montageNr].secondContactGlobalIdx, lStartSampleRead, EDFSEEK_SET);//set start Sample
			_int64 outA = edfread_physical_samples((int)edfhdr.handle, (int)m_montages[montageNr].firstContactGlobalIdx, (int)lNumberSamplesRead, bufA);
			_int64 outB = edfread_physical_samples((int)edfhdr.handle, (int)m_montages[montageNr].secondContactGlobalIdx, (int)lNumberSamplesRead, bufB);

			if (outA == -1 || outB == -1)
				return false;

			for (long i = 0; i < lNumberSamplesRead; ++i) {
				reformattedSamples[montageNr][i] = bufA[i] - bufB[i];
			}
			delete[] bufA;
			delete[] bufB;

			edfclose_file(edfhdr.handle);
		}
	}

	edf_hdr_struct edfhdr;
	_int64 headerReadOK = edfopen_file_readonly(getFileName().c_str(), &edfhdr, EDFLIB_DO_NOT_READ_ANNOTATIONS);
	edfclose_file(edfhdr.handle);
	if (headerReadOK == -1)
		return false;


	std::string writePath = m_patientPath +  m_patientName + "_mtg.edf";
	remove(writePath.c_str());
	int writeHandle = edfopen_file_writeonly(writePath.c_str(), EDFLIB_FILETYPE_EDFPLUS, m_montages.size());

	//Write header data specific to file
	double dataRecordDurationInSec = (edfhdr.datarecord_duration * 100.0) / (1000.0 * 1000.0 * 1000.0);
	edf_set_startdatetime( writeHandle, edfhdr.startdate_year, edfhdr.startdate_month, edfhdr.startdate_day, edfhdr.starttime_hour, edfhdr.starttime_minute, edfhdr.starttime_second);
	edf_set_patientname( writeHandle, m_patientName.c_str());
	edf_set_patientcode( writeHandle, m_patientName.c_str());
	edf_set_patient_additional( writeHandle, edfhdr.patient_additional);
	edf_set_admincode( writeHandle, edfhdr.admincode);
	edf_set_technician( writeHandle, edfhdr.technician);
	edf_set_equipment( writeHandle, edfhdr.equipment);
	edf_set_recording_additional( writeHandle, edfhdr.recording_additional);
	for (int nrEDF_Signal = 0; nrEDF_Signal < m_montages.size(); ++nrEDF_Signal) {
		//Write header data specific to each signal
		double samplingRate = (double)edfhdr.signalparam[nrEDF_Signal].smp_in_datarecord / dataRecordDurationInSec;
		edf_set_samplefrequency( writeHandle,  nrEDF_Signal, samplingRate);
		edf_set_datarecord_duration(writeHandle, 100000);
		edf_set_physical_maximum( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].phys_max);
		edf_set_physical_minimum( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].phys_min);
		edf_set_digital_maximum( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].dig_max);
		edf_set_digital_minimum( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].dig_min);
		edf_set_label( writeHandle,  nrEDF_Signal, m_montages[nrEDF_Signal].montageName.c_str());
		edf_set_prefilter( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].prefilter);
		edf_set_transducer( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].transducer);
		edf_set_physical_dimension( writeHandle,  nrEDF_Signal, edfhdr.signalparam[nrEDF_Signal].physdimension);
	}

	int blockIdx = 1;
	int blockSize = (double)edfhdr.signalparam[0].smp_in_datarecord / dataRecordDurationInSec;
	while ((blockIdx * blockSize) < lNumberSamplesRead) {

		//Write signals
		for (int nrEDF_Signal = 0; nrEDF_Signal < m_montages.size(); ++nrEDF_Signal) {
			double samplingRate = (double)edfhdr.signalparam[nrEDF_Signal].smp_in_datarecord / dataRecordDurationInSec;
			double *writeBuf = new double[samplingRate];

			for (int i = 0; i < samplingRate; ++i) {
				long long globalSmplIdx = (blockIdx - 1) * blockSize + i;
				writeBuf[i]= reformattedSamples[nrEDF_Signal][globalSmplIdx];
			}

			edfwrite_physical_samples(writeHandle, writeBuf);
			delete[] writeBuf;
		}

		blockIdx++;
	}
	

	edfclose_file(writeHandle);


	return true;
}

std::string IO_EDFPlus_SignalImporter_c::getPatientNameAndPath(void) {
	for (int i = m_strFilename.length() - 1; i >= 0; i--) {
		if (m_strFilename[i] == 92) {
			m_patientName = m_strFilename.substr(i + 1, m_strFilename.length() - i);
			m_patientPath = m_strFilename.substr(0, i + 1);
			break;
		}
	}

	m_patientName.pop_back(); m_patientName.pop_back(); m_patientName.pop_back(); m_patientName.pop_back();

	return m_patientName;
}