#include <string>
#include <vector>
#include <fstream>
#include <direct.h>

#include "Plotter.h"

//#define DEFINE_PLOTTER

#ifdef DEFINE_PLOTTER
#include "T:\VisualStudio\OscillationsDetector\MOSDET\gnuplot-iostream-master\gnuplot-iostream.h"
#endif

void Plotter::plotXY(std::vector<double> &x, std::vector<double> &y, std::string &wdwTitle, double xtics, double ytics) {

#ifdef DEFINE_PLOTTER

	Gnuplot plot("gnuplot -persist");

	std::string imageFilename = ( wdwTitle + ".png");

	std::vector<std::pair<double, double> > xy_pts_A;
	for (int i = 0; i < x.size(); ++i) {
		xy_pts_A.push_back(std::make_pair(x.at(i), y.at(i)));
	}
	plot << "reset\n";
	//plot << "set term wxt " << m_windowIdx << " title '" << wdwTitle << "'\n";
	//plot << "set term svg size 1920, 1200\n";
	plot << "set terminal pngcairo size 1920, 1200 enhanced font 'Verdana,10'\n";
	plot << "set output\"" << imageFilename << "\"\n";
	plot << "set title '" << wdwTitle << "' font \"arial, 20\"\n";
	plot << "set autoscale y\n";
	plot << "set autoscale x\n";
	plot << "set xrange[" << std::to_string(*std::min_element(x.begin(), x.end())) << ":" << std::to_string(*std::max_element(x.begin(), x.end())) << "]\n";
	//plot << "set xrange[0:30]\n";
	//plot << "set format y \"%.1f\"\n";
	plot << "set yrange[" << std::to_string(2*(*std::min_element(y.begin(), y.end()))) << ":" << std::to_string(2 * (*std::max_element(y.begin(), y.end()))) << "]\n";
	plot << "set key box opaque\n";
	plot << "set ylabel 'Amplitude (V)'\n";
	plot << "set xlabel 'Samples'\n";
	plot << "set grid lw 0.3, lw 0.3\n";
	plot << "set tics\n";
	plot << "set grid ytics\n";
	if (xtics != 0)
		plot << "set xtics " << xtics << "\n";
	if (ytics != 0)
		plot << "set ytics " << ytics << "\n";
	//plot << "set mxtics 1\n";
	plot << "plot '-' with lines\n";
	plot.send1d(xy_pts_A);
	plot << "unset ylabel\n";
	plot << "unset ytics\n";
	plot << "unset xlabel\n";

#else
	writeVectorsToFile(x, y, wdwTitle);
#endif


}

bool Plotter::Plot_3D_Profile(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::string &wdwTitle, double xtics, double xmin, double xmax, double ytics, double ymin, double ymax) {

	/*unsigned dataSize = x.size();
	std::string filename = ("FilesOutput/" + wdwTitle + ".dat");
	const char * filenameC = filename.c_str();
	remove(filenameC);
	std::ofstream data;
	data.open(filenameC);
	data << "#X Y Z\n";
	for (unsigned i = 0; i < dataSize; ++i) {
		data << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
		if ((i + 1 < dataSize) && (x[i] != x[i + 1]))
			data << "\n";
	}
	data << "\n";
	data.close();//close file

	//return true;

	Gnuplot plot("gnuplot -persist");
	std::string imageFilename = ("FilesOutput/" + wdwTitle + ".png");
	plot << "reset\n";
	//plot << "set terminal wxt size 1200, 800 enhanced font 'Verdana,10' persist\n";
	plot << "set terminal pngcairo size 1920, 1200 enhanced font 'Verdana,10'\n";
	plot << "set output\"" << imageFilename << "\"\n";
	//plot << "unset key\n";

	//border
	plot << "set style line 11 lc rgb '#808080' lt 1\n";
	plot << "set border 3 front ls 11\n";
	plot << "set tics nomirror out scale 0.75\n";
	if (xtics != 0)
		plot << "set xtics " << xtics << "\n";
	if (ytics != 0)
		plot << "set ytics " << ytics << "\n";

	//Colorbar
	//disable colorbar tics
	plot << "set cbtics scale 0\n";
	bool useDefaultPalette = true;
	if (useDefaultPalette) {
		//Default matlab palette colors
		plot << "set palette defined(0 \"#000090\","
			<< "1 \"#000fff\","
			<< "2 \"#0090ff\","
			<< "3 \"#0fffee\","
			<< "4 \"#90ff70\","
			<< "5 \"#ffee00\","
			<< "6 \"#ff7000\","
			<< "7 \"#ee0000\","
			<< "8 \"#7f0000\")"
			<< "\n";
	}
	else {
		//New Matlab palette colors
		std::string matlabColormap = ("gnuplot-iostream-master/parula.pal");
		plot << "load '" << matlabColormap << "'\n";
	}


	plot << "set pm3d map\n";
	plot << "set dgrid3d 200, 200\n";
	//plot << "set hidden3d\n";
	if (xmin != 0 || xmax != 0)
		plot << "set xrange[" << xmin << ":" << xmax << "]\n";
	if (ymin != 0 || ymax != 0)
		plot << "set yrange[" << ymin << ":" << ymax << "]\n";


	plot << "splot \"" << filename << "\" using 1:2:3\n";*/

	return true;
}

void Plotter::writeVectorsToFile(std::vector<double> &x, std::vector<double> &y, std::string &filename) {
	filename += ".txt";
	long nrSamples = x.size();
	remove(filename.c_str());

	std::ofstream vectorsFile;
	vectorsFile.open(filename, std::ios::app);    // open file for appending
	vectorsFile << "VectorA\t" << "VectorB\n";

	for (long i = 0; i < nrSamples; i++) {
		vectorsFile << x[i] << "\t" << y[i] << "\n";
	}

	vectorsFile.close();
}