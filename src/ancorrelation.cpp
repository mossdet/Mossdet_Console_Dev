// Based on the book numerical recipes, The Art of Scientific Computing, Third Edition

#include <math.h>
#include <limits>
#include <assert.h>
#include "ancorrelation.h"
#include "HeapSorter.h"


spearRslt Correlation_c::getSpearmanCorrCoeff(std::vector<double> &A, std::vector<double> &B, long& ccLag) {

	// Recalculate vector ranks for each lag
	spearRslt result;
	std::vector<double> signalA;
	std::vector<double> signalB;


	/*if (ccLag >= 0) { //Positive lag moves start point of A to Asample number lag (equvalent to rotating A to the left with B static)
		for (unsigned i = 0; i < B.size() - ccLag; ++i)
			signalB.push_back(B[i]);
		for (unsigned i = ccLag; i < A.size(); ++i)
			signalA.push_back(A[i]);
	} else {//Negative lag moves start point of B to Bsample number lag (equvalent to rotating B to the left with A static)
		for (unsigned i = 0; i < A.size() - (-1 * ccLag); ++i)
			signalA.push_back(A[i]);
		for (unsigned i = (-1 * ccLag); i < B.size(); ++i)
			signalB.push_back(B[i]);
	}*/

	unsigned sampleNr = A.size();
	//Negative lag moves start point of A to Asample number lag (equvalent to rotating A to the left with B static)
	//Positive lag moves start point of B to Bsample number lag (equvalent to rotating B to the left with A static)

	for (int i = 0; i < sampleNr; i++) {
		if (i + ccLag < 0)
			continue;
		if (i + ccLag >= sampleNr)
			break;
		signalB.push_back(B[i + ccLag]);
		signalA.push_back(A[i]);
	}

	assert(signalA.size() == signalB.size());

	spear(signalA, signalB, result.d, result.zd, result.probd, result.rs, result.probrs);

	return result;

	// Do not recalculate vector ranks for each lag
	/*spearRslt result;
	std::vector<double> signalA = A;
	std::vector<double> signalB = B;

	double sA, sB;
	double aAvg = 0, bAvg = 0;
	double aSqrSum = 0, bSqrSum = 0;
	double num = 0, denom = 1;
	int sampleNr = signalA.size();


	// Get Ranks
	Sorter::heapSort2(signalA, signalB);
	crank(signalA, sA);

	Sorter::heapSort2(signalB, signalA);
	crank(signalB, sB);

	result.d = result.zd = result.probd = result.probrs = 0;
	result.rs = getCorrCoeff(signalA, signalB, ccLag);

	return result;*/
}

double Correlation_c::getCorrCoeff(std::vector<double>& signalA, std::vector<double>& signalB, long& ccLag) {

	double aAvg = 0, bAvg = 0;
	double aSqrSum = 0, bSqrSum = 0;
	double num = 0, denom = 1;
	int sampleNr = signalA.size();

	for (int i = 0; i < sampleNr; i++) {
		aAvg += signalA[i];
		bAvg += signalB[i];
	}
	aAvg /= sampleNr;
	bAvg /= sampleNr;


	for (int i = 0; i < sampleNr; i++) {
		if (i + ccLag < 0)
			continue;
		if (i + ccLag >= sampleNr)
			break;
		num += (signalA[i] - aAvg)*(signalB[i + ccLag] - bAvg);
		aSqrSum += pow(signalA[i] - aAvg, 2);
		bSqrSum += pow(signalB[i + ccLag] - bAvg, 2);
	}
	denom = sqrt(aSqrSum*bSqrSum);
	if (denom == 0) {
		return 0;
	}
	return num / denom;
}

void Correlation_c::crank(std::vector<double> &w, double &s) { // Given a sorted array w[1..n], replaces the elements by their rank, including midranking of ties, and returns as s the sum of f3 ? f, w here f is the number of elements in each tie.
	unsigned n = w.size();
	double t, rank = 1;
	s = 0.0;
	unsigned idx = 0, idxTwo = 0, idxThree = 0;
	while (idx < n - 1) {
		if (w[idx + 1] != w[idx]) { //Not a tie.
			w[idx] = rank;
			rank++;
			idx++;
		}
		else { //A tie :
			for (idxTwo = idx + 1; idxTwo < n && w[idxTwo] == w[idx]; idxTwo++);		// How far does it go ?

			rank = ((idxTwo - 1 + 1) + (idx + 1)) / 2;									// This is the mean rank of the tie, so enter it into all the tied entries, and update s.
			for (idxThree = idx; idxThree < idxTwo; idxThree++)
				w[idxThree] = rank;
			t = idxTwo - idx;
			s += t*t*t - t;
			idx = idxTwo;
		}
	}
	if (idx == n - 1)
		w[n - 1] = n;	//If the last element was not tied, this is its rank.
}

double Correlation_c::erfcc(const double x) //Returns the complementary error function erfc.x / with fractional error everywhere less than 1.2 x 10exp-7.
{
	double t, z = fabs(x), ans;
	t = 2. / (2. + z);
	ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
		t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
			t*(-0.82215223 + t*0.17087277)))))))));
	return (x >= 0.0 ? ans : 2.0 - ans);
}

void Correlation_c::swapElements(double &a, double &b) {
	double swap = a;
	a = b;
	b = swap;
}


void Correlation_c::pearsn(std::vector<double> &x, std::vector<double> &y, double &r, double &prob, double &z) {
	const double TINY = 1.0e-20;													// Will regularize the unusual case of complete correlation.
	Beta beta;
	int j, n = x.size();
	double yt, xt, t, df;
	double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;
	for (j = 0; j<n; j++) {															//Find the means.
		ax += x[j];
		ay += y[j];
	}
	ax /= n;
	ay /= n;
	for (j = 0; j<n; j++) {															//Compute the correlation coefficient.
		xt = x[j] - ax;
		yt = y[j] - ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	r = sxy / (sqrt(sxx*syy) + TINY);
	z = 0.5*log((1.0 + r + TINY) / (1.0 - r + TINY));									//Fisher’s z transformation.
	df = n - 2;
	t = r*sqrt(df / ((1.0 - r + TINY)*(1.0 + r + TINY)));								//Equation(14.5.5).
	prob = beta.betai(0.5*df, 0.5, df / (df + t*t));									//Student’s t probability.
	prob = erfcc(abs(z*sqrt(n - 1.0)) / 1.4142136);
	//For large n, this easier computation of prob, using the short routine erfcc, would give approximately the same value.
}



/*Given two data arrays, data1[1..n] and data2[1..n], this routine returns their sum - squared
difference of ranks as D, the number of standard deviations by which D deviates from its nullhypothesis
expected value as zd, the two - sided significance level of this deviation as probd,
Spearman’s rank correlation rs as rs, and the two - sided significance level of its deviation from
zero as probrs.The external routines crank(below) and sort2(§8.2) are used.A small value
of either probd or probrs indicates a significant correlation(rs positive) or anticorrelation
(rs negative).*/
void Correlation_c::spear(std::vector<double> &a, std::vector<double> &b, double &d, double &zd, double &probd, double &rs, double &probrs) {
	assert(a.size() == b.size());
	assert(a.size() >0 && b.size() > 0);
	Beta beta;
	double vard, t, sg, sf, fac, en3n, en, df, aved;

	std::vector<double> wksp1 = a;
	std::vector<double> wksp2 = b;

	//Sort each of the data arrays, and convert the entries to ranks.The values sf and sg return the sums SUM(fk^3 - fk) and SUM(gm^3 - gm), respectively.
	Sorter::heapSort2(wksp1, wksp2);
	crank(wksp1, sf);

	//Sort each of the data arrays, and convert the entries to ranks.The values sf and sg return the sums SUM(fk^3 - fk) and SUM(gm^3 - gm), respectively.
	Sorter::heapSort2(wksp2, wksp1);
	crank(wksp2, sg);

	d = 0.0;
	for (unsigned j = 0; j < wksp1.size(); j++) //Sum the squared difference of ranks.
		d += SQR(wksp1[j] - wksp2[j]);
	en = wksp1.size();
	en3n = en*en*en - en;
	aved = en3n / 6.0 - (sf + sg) / 12.0; //Expectation value of D,
	fac = (1.0 - sf / en3n)*(1.0 - sg / en3n);
	vard = ((en - 1.0)*en*en*SQR(en + 1.0) / 36.0)*fac; //and variance of D give
	zd = (d - aved) / sqrt(vard); //number of standard deviations and significance.
	probd = erfcc(fabs(zd) / 1.4142136);
	rs = (1.0 - (6.0 / en3n)*(d + (sf + sg) / 12.0)) / sqrt(fac); //Rank correlation coefficient,
	fac = (rs + 1.0)*(1.0 - (rs));
	if (fac > 0.0) {
		t = (rs)*sqrt((en - 2.0) / fac); //and its t value,
		df = en - 2.0;
		probrs = beta.betai(0.5*df, 0.5, df / (df + t*t)); //give its significance.
	}
	else
		probrs = 0.0;
}


double Beta::betai(const double a, const double b, const double x) {
	//Returns incomplete beta function Ix.a; b / for positive a and b, and x between 0 and 1.
	double bt;
	if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
	if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) return x;
	if (a > SWITCH && b > SWITCH) return betaiapprox(a, b, x);
	bt = exp(gm.gammln(a + b) - gm.gammln(a) - gm.gammln(b) + a*log(x) + b*log(1.0 - x));
	if (x < (a + 1.0) / (a + b + 2.0)) return bt*betacf(a, b, x) / a;
	else return 1.0 - bt*betacf(b, a, 1.0 - x) / b;
}

double Beta::betacf(const double a, const double b, const double x) { //Evaluates continued fraction for incomplete beta function by modified Lentz’s method (5.2).User should not call directly.
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;
	qab = a + b;			//These q’s will be used in factors that
	qap = a + 1.0;			//occur in the coefficients(6.4.6).
	qam = a - 1.0;
	c = 1.0;				//First step of Lentz’s method.
	d = 1.0 - qab*x / qap;
	if (fabs(d) < FPMIN) d = FPMIN;
	d = 1.0 / d;
	h = d;
	for (m = 1; m<10000; m++) {
		m2 = 2 * m;
		aa = m*(b - m)*x / ((qam + m2)*(a + m2));
		d = 1.0 + aa*d;								//One step(the even one) of the recurrence.
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		h *= d*c;
		aa = -(a + m)*(qab + m)*x / ((a + m2)*(qap + m2));
		d = 1.0 + aa*d;										//Next step of the recurrence(the odd one)
		if (fabs(d) < FPMIN) d = FPMIN;
		c = 1.0 + aa / c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0 / d;
		del = d*c;
		h *= del;
		if (fabs(del - 1.0) <= EPS) break;					//Are we done ?
	}
	return h;
}

double Beta::betaiapprox(double a, double b, double x) {
	//Incomplete beta by quadrature.Returns Ix.a; b / .User should not call directly.
	int j;
	double xu, t, sum, ans;
	double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
	double lnmu = log(mu), lnmuc = log(1. - mu);
	t = sqrt(a*b / (SQR(a + b)*(a + b + 1.0)));
	if (x > a / (a + b)) {
		//Set how far to integrate into the tail :
		if (x >= 1.0) return 1.0;
		xu = MIN(1., MAX(mu + 10.*t, x + 5.0*t));
	}
	else {
		if (x <= 0.0) return 0.0;
		xu = MAX(0., MIN(mu - 10.*t, x - 5.0*t));
	}
	sum = 0;
	for (j = 0; j<18; j++) {
		//Gauss - Legendre.
		t = x + (xu - x)*Gauleg18::y[j];
		sum += Gauleg18::w[j] * exp(a1*(log(t) - lnmu) + b1*(log(1 - t) - lnmuc));
	}
	ans = sum*(xu - x)*exp(a1*lnmu - gm.gammln(a) + b1*lnmuc - gm.gammln(b) + gm.gammln(a + b));
	return ans>0.0 ? 1.0 - ans : -ans;
}

double Beta::invbetai(double p, double a, double b) {
	//Inverse of incomplete beta function. Returns x such that Ix.a; b/ D p for argument p between 0 and 1.
	const double EPS = 1.e-8;
	double pp, t, u, err, x, al, h, w, afac, a1 = a - 1., b1 = b - 1.;
	int j;
	if (p <= 0.) return 0.;
	else if (p >= 1.) return 1.;
	else if (a >= 1. && b >= 1.) {								//Set initial guess. See text.
		pp = (p < 0.5) ? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753 + t*0.27061) / (1. + t*(0.99229 + t*0.04481)) - t;
		if (p < 0.5) x = -x;
		al = (SQR(x) - 3.) / 6.;
		h = 2. / (1. / (2.*a - 1.) + 1. / (2.*b - 1.));
		w = (x*sqrt(al + h) / h) - (1. / (2.*b - 1) - 1. / (2.*a - 1.))*(al + 5. / 6. - 2. / (3.*h));
		x = a / (a + b*exp(2.*w));
	}
	else {
		double lna = log(a / (a + b)), lnb = log(b / (a + b));
		t = exp(a*lna) / a;
		u = exp(b*lnb) / b;
		w = t + u;
		if (p < t / w) x = pow(a*w*p, 1. / a);
		else x = 1. - pow(b*w*(1. - p), 1. / b);
	}
	afac = -gm.gammln(a) - gm.gammln(b) + gm.gammln(a + b);
	for (j = 0; j<10; j++) {
		if (x == 0. || x == 1.) return x;						//a or b too small for accurate calculation.
		err = betai(a, b, x) - p;
		t = exp(a1*log(x) + b1*log(1. - x) + afac);
		u = err / t; Halley:
		x -= (t = u / (1. - 0.5*MIN(1., u*(a1 / x - b1 / (1. - x)))));
		if (x <= 0.) x = 0.5*(x + t);							//Bisect if x tries to go neg or > 1.
		if (x >= 1.) x = 0.5*(x + t + 1.);
		if (fabs(t) < EPS*x && j > 0) break;
	}
	std::numeric_limits<double>::min;
	return x;
}