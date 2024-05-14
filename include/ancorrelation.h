
#pragma once

#include <vector>

typedef struct {
	double d;
	double zd;
	double probd;
	double rs;
	double probrs;
} spearRslt;

static float sqrarg;
//static double swap;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 :sqrarg*sqrarg)
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
//#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
//#define FPMIN 1.0e-30				// Close to smallest representable floating-point number.
//#define EPS 6.0e-8				// Relative error, or absolute error near a zero of Ci(x).

class Correlation_c {
public:
	static spearRslt getSpearmanCorrCoeff(std::vector<double> &signalA, std::vector<double> &signalB, long& ccLag);
	static double getCorrCoeff(std::vector<double>& signalA, std::vector<double>& signalB, long& ccLag);
	static void pearsn(std::vector<double> &x, std::vector<double> &y, double &r, double &prob, double &z);
	static void spear(std::vector<double> &a, std::vector<double> &b, double &d, double &zd, double &probd, double &rs, double &probrs);
	static void crank(std::vector<double> &w, double &s);
	static double erfcc(const double x);
	static void swapElements(double &a, double &b);

};

class Gauleg18 {
public:
	//Abscissas and weights for Gauss - Legendre quadrature.
	static const int ngau = 18;
	//static const double y[18];
	//static const double w[18];

	const double y[18] = { 0.0021695375159141994,
		0.011413521097787704,0.027972308950302116,0.051727015600492421,
		0.082502225484340941, 0.12007019910960293,0.16415283300752470,
		0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
		0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
		0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
		0.87126389619061517, 0.95698180152629142 };

	const double w[18] = { 0.0055657196642445571,
		0.012915947284065419,0.020181515297735382,0.027298621498568734,
		0.034213810770299537,0.040875750923643261,0.047235083490265582,
		0.053244713977759692,0.058860144245324798,0.064039797355015485,
		0.068745323835736408,0.072941885005653087,0.076598410645870640,
		0.079687828912071670,0.082187266704339706,0.084078218979661945,
		0.085346685739338721,0.085983275670394821 };
};

class Gamma {
public:
	double gammln(const double xx) {
		//Returns the value lnOE.xx / for xx > 0.
		int j;
		double x, tmp, y, ser;
		static const double cof[14] = { 57.1562356658629235,-59.5979603554754912,
			14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
			.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
			-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
			.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5 };
		if (xx <= 0) throw("bad arg in gammln");
		y = x = xx;
		tmp = x + 5.24218750000000000;				//Rational 671 / 128.
		tmp = (x + 0.5)*log(tmp) - tmp;
		ser = 0.999999999999997092;
		for (j = 0; j<14; j++) ser += cof[j] / ++y;
		return tmp + log(2.5066282746310005*ser / x);
	}
};

class Beta : Gauleg18 {
public:

	//Object for incomplete beta function.Gauleg18 provides coefficients for Gauss - Legendre quadrature.
	static const int SWITCH = 3000; //When to switch to quadrature method.
	const double EPS = 2.2204460492503131e-16; // std::numeric_limits<double>::epsilon();
	const double FPMIN = 1.0020841800044864e-292;//1.0e-30;//std::numeric_limits<double>::min() / Beta::EPS;

	Gamma gm;
	double betai(const double a, const double b, const double x);

	double betacf(const double a, const double b, const double x);

	double betaiapprox(double a, double b, double x);

	double invbetai(double p, double a, double b);
};

class Betadist : Beta {
	//Beta distribution, derived from the beta function Beta.
	double alph, bet, fac;
	Betadist(double aalph, double bbet) : alph(aalph), bet(bbet) {
		//Constructor.Initialize with alpha and beta
		if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Betadist");
		fac = gm.gammln(alph + bet) - gm.gammln(alph) - gm.gammln(bet);
	}
	double p(double x) {
		//Return probability density function.
		if (x <= 0. || x >= 1.) throw("bad x in Betadist");
		return exp((alph - 1.)*log(x) + (bet - 1.)*log(1. - x) + fac);
	}
	double cdf(double x) {
		//Return cumulative distribution function.
		if (x < 0. || x > 1.) throw("bad x in Betadist");
		return betai(alph, bet, x);
	}
	double invcdf(double p) {
		//Return inverse cumulative distribution function.
		if (p < 0. || p > 1.) throw("bad p in Betadist");
		return invbetai(p, alph, bet);
	}
};
