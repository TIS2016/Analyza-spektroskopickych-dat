#pragma once

#define _USE_MATH_DEFINES // for C++ pi=M_PI,       1/pi=M_1_PI

/* Definitions of useful mathematical constants
* M_E        - e
* M_LOG2E    - log2(e)
* M_LOG10E   - log10(e)
* M_LN2      - ln(2)
* M_LN10     - ln(10)
* M_PI       - pi
* M_PI_2     - pi/2
* M_PI_4     - pi/4
* M_1_PI     - 1/pi
* M_2_PI     - 2/pi
* M_2_SQRTPI - 2/sqrt(pi)
* M_SQRT2    - sqrt(2)
* M_SQRT1_2  - 1/sqrt(2)
*/

//Pi
#define M_1oSqrt2PI  (1 / sqrt(2 * M_PI))
#define M_1oSqrtPI  (1 / sqrt(M_PI))
#define M_SqrtPI  (sqrt(M_PI))
#define M_2SqrtPI  (2 * sqrt(M_PI))

//Ln2...
#define M_SqrtLn2 (sqrt(M_LN2))
#define M_Sqrt2Ln2 (sqrt(2* M_LN2))
#define M_1o2Ln2 (1 / (2* M_LN2))
#define M_1oSqrt2Ln2 (1 / sqrt(2* M_LN2))

//Pi&Ln
#define M_SqrtLn2oPi (sqrt(M_LN2/M_PI))

//For Voigt calc.
#define CVOIGT 1.12837916709551257

//For other stuff

#define P_C 299792458
#define P_1o2C (1/(2.0*299792458))
#define P_CSqrtLn2 (299792458 * sqrt(M_LN2))
#define P_CSqrtPi (299792458 * sqrt(M_PI))
#define P_CoLn2 (P_C/M_LN2)
#define P_CoSqrtLn2 (P_C/sqrt(M_LN2))

#include <cmath>
#include <complex>
#include <vector>
#include <ppl.h>
#include "concurrent_vector.h"

using namespace concurrency;

// Pour les nombres de complexes
typedef std::complex<double> dcplx;


/***************************************************************************/
//Main Function
double PeakAdv(double xi, double *Par, bool am);

//Peak Profiles
double LorentzA(const double xi, const double x0, const double a, const double w0, const double d0);
double LorentzS(const double xi, const double x0, const double s, const double w0, const double d0);

double GaussA(const double xi, const double x0, const double a, const double wD, const double d0);
double GaussS(const double xi, const double x0, const double s, const double wD, const double d0);

double VoigtHA(const double xi, const double x0, const double a, const double wD, const double w0, const double d0);
double VoigtHS(const double xi, const double x0, const double s, const double wD, const double w0, const double d0);

double VoigtA(const double xi, const double x0, const double a, const double wD, const double w0, const double d0);
double VoigtS(const double xi, const double x0, const double s, const double wD, const double w0, const double d0);

double HTPA(double &xi, double *Par);

//To calc. the complex err. function w(z)

void cwerf(double x, double y, double *vr, double *vi);

dcplx complexErrFunc(const dcplx &z);
dcplx CPF(const dcplx &z);
dcplx CPF3(const dcplx &z);


class simulation {

public:

	/***********************/
	/* Templates de calcul */
	/***********************/

	template <class T>
	inline bool isnan(T s)
	{
		// By IEEE 754 rule, 2*Inf equals Inf
		return (s != s);
	}

	// Calcul d'une valeur inverse
	template <class T>
	inline T fcnInverse(T x) const
	{
		return ((T)1.0 / x);
	};
};