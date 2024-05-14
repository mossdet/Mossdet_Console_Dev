#include <math.h>
#include <fstream>
#include <string>

#include "WaveletTransform.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define sqrx(x) ((x<=1e-36)?0:sqrt(x)) // safe square root

WaveletTransform::WaveletTransform(double samplingRate, double startFreq, double endFreq, double deltaFreq, int nrOscills)
{
	m_dSamplingRate = samplingRate;
	m_dStartFrequency = startFreq; 40;
	m_dEndFrequency = endFreq; 560;
	m_dDeltaFrequency = deltaFreq; 40;
	m_nrOscillations = nrOscills; 4;

}


WaveletTransform::~WaveletTransform()
{
}

bool WaveletTransform::run(matrixStd &data) {
	// get number components
	unsigned int nComponents = getNumberComponents();
	nComponents = (m_dEndFrequency - m_dStartFrequency) / m_dDeltaFrequency;
	if (nComponents == 0){
		return false;
	}
	long nrChanns = data.size();
	long nrParamRows;
	long nrColumns;

	// allocate result data
	nrParamRows = data.size()* (m_dEndFrequency - m_dStartFrequency) / m_dDeltaFrequency;
	nrColumns = data[0].size();

	m_Parameters.resize(nrParamRows);
	for (long r = 0; r < nrParamRows; ++r)
		m_Parameters[r].resize(nrColumns);

	// loop over channels
	for (int h = 0; h < nrChanns; h++) {
		// loop over wavelet components
		for (long nc = 0; nc < nComponents; nc++){
			// get convolution kernel
			double dCenterFrequency = m_dStartFrequency + nc * m_dDeltaFrequency;
			std::vector<double>  kernel = getConvolutionKernel(dCenterFrequency);

			// loop over samples
#pragma omp parallel for// private (data, kernel)
			for (long l = 0; l < nrColumns; l++) {
				// perform convolution
				long lKernelLength = kernel.size();

				if (l - lKernelLength < 0)
					lKernelLength = l;

				for (long m = 0; m < lKernelLength; m++) {
					m_Parameters[h  * nComponents + nc][l] = m_Parameters[h * nComponents + nc][l] + kernel[m] * data[h][l - m];
				}
			}
		}

	}

	return true;
}

unsigned int WaveletTransform::getNumberComponents() const{
	unsigned int nNumberComponents = 0;
	// preform checks
	if (m_dDeltaFrequency < 1E-20 || m_dDeltaFrequency >(m_dEndFrequency - m_dStartFrequency))
		return nNumberComponents;

	nNumberComponents = (m_dEndFrequency - m_dStartFrequency) / m_dDeltaFrequency;

	return nNumberComponents;
}

std::vector<double> WaveletTransform::getConvolutionKernel(double dCenterFrequency){
	std::vector<double> kernel;
	// check for valid sampling frequency
	if (m_dSamplingRate < 1e-12 || dCenterFrequency < 1e-12)
		return kernel;

	// compute relative frequency for sinus function
	double dRelFrequency = dCenterFrequency / m_dSamplingRate;

	unsigned int nrCycles = m_nrOscillations;

	// determine kernel length depending on sampling frequency and center frequency
	long lkernelLenght = (nrCycles / dCenterFrequency) * m_dSamplingRate;
	kernel.resize(lkernelLenght);

	double d2PI = 2 * M_PI;
	double dGausScale = 1 / dRelFrequency * nrCycles / d2PI;
	double dGausScaleSqr = dGausScale * dGausScale;
	double dAmpScale = 1.0 / sqrx(1 / dRelFrequency * nrCycles / d2PI)*1.0 / sqrx(sqrx(M_PI));
	for (long k = 0; k < lkernelLenght; k++){
		kernel[k] = cos(dRelFrequency * k * d2PI)* dAmpScale *exp(-1.0* double((k - lkernelLenght / 2)*(k - lkernelLenght / 2) / dGausScaleSqr) / 2.0);
	}

	/*std::string name = "Kernel" + std::to_string((int)dCenterFrequency) + "Hz-CenterFrequency-" + std::to_string(nrCycles) + "Cycles";
	std::string filename = "T:\\VisualStudio\\OscillationsDetector\\MOSDET\\" + name + ".dat";
	std::ofstream kernelDataFile;
	kernelDataFile.open(filename);
	kernelDataFile << "Kernel\n";
	for (long i = 0; i < kernel.size(); i++)
		kernelDataFile << kernel [i] << "\n";
	kernelDataFile.close();*/

	return kernel;
}

bool WaveletTransform::getResult(matrixStd &outParameters) {
	if (m_Parameters.size() == 0)
		return false;
	outParameters = m_Parameters;
	return true;
}