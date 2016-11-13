/*Peter selection of LinesM and Serge Simulation class - ready for parallel computing*/

#include "linesAdv.h"
#include "Cerf.h"

/*********Main Function****************************************/

double PeakAdv(double xi, double *Par, bool am)
{
	//read peak shape
	int shape = *Par;

	switch (shape)
	{
		//lorentz***********************
	case 0:
		if (am) return LorentzA(xi, *(Par + 1), *(Par + 2), *(Par + 5), *(Par + 6));
		else return LorentzS(xi, *(Par + 1), *(Par + 3), *(Par + 5), *(Par + 6));

		break;

		//gauss***********************
	case 1:
		if (am) return GaussA(xi, *(Par + 1), *(Par + 2), *(Par + 4), *(Par + 6));
		else return GaussS(xi, *(Par + 1), *(Par + 3), *(Par + 4), *(Par + 6));

		break;

		//voigt***********************
	case 2:
		if (am) return VoigtA(xi, *(Par + 1), *(Par + 2), *(Par + 4), *(Par + 5), *(Par + 6));
		else return VoigtS(xi, *(Par + 1), *(Par + 3), *(Par + 4), *(Par + 5), *(Par + 6));

		break;

		//voigtH***********************
	case 3:
		if (am) return VoigtHA(xi, *(Par + 1), *(Par + 2), *(Par + 4), *(Par + 5), *(Par + 6));
		else return VoigtHS(xi, *(Par + 1), *(Par + 3), *(Par + 4), *(Par + 5), *(Par + 6));

		break;

		//HTP***********************
	case 8:
		if (am) return HTPA(xi, Par);
		else return HTPA(xi, Par);
	
	default:
		break;
	}
	return 0.0;
}

/*in Par array the parameters are stored in following order
Peak param:
0 shape
1 center
2 amp.
3 surf.
4 wD            HWHM
5 w0            Collisional HWHM           
6 d0            defined as x0+d0
7 w2
8 d2
9 eta
10 nuVC
11 beta			for Galatry
12 vMP
13 sample
14 center_in
...
in total = ParNumber
*/

/*Part 1: Concurence stuff*/

/*
Lorentzian->lorentz() A or S according to the amplitude/surface leading parameter

computed as: 
where dx = x - x0, so g = gamma and g = HWHM, a- amplitude, s-surface

a g^2 / (dx^2 + g^2), 
s/pi * g / (dx^2 + g^2)

*/

double LorentzA(const double xi, const double x0, const double a, const double w0, const double d0)
{
	return a * (w0 * w0)/ ((w0 * w0) + (xi-(x0 + d0))*(xi - (x0 + d0)));
}

double LorentzS(const double xi, const double x0, const double s, const double w0, const double d0)
{
	return (s * M_1_PI) * w0 / ((w0 * w0) + (xi - (x0 + d0))*(xi - (x0 + d0)));
}

/*
Gaussian->gauss() A or S...
computed as

... wD = HWHM

*/

double GaussA(const double xi, const double x0, const double a, const double wD, const double d0)
{
	return a * exp(-(xi - (x0 + d0))*(xi - (x0 + d0))/(2 * M_1o2Ln2 * wD * wD));
}

double GaussS(const double xi, const double x0, const double s, const double wD, const double d0)
{
	return s * (M_SqrtLn2oPi / wD) * exp(-(xi - (x0 + d0))*(xi - (x0 + d0)) / (2 * M_1o2Ln2 * wD * wD));
}

/*
Voigt peak

Vs(x; wD,w0) = SQRT(ln2/PI) * Re[w(z)]/wD
z = SQRT(ln2) * (x/wD) , SQRT(ln2) * (w0/wD)

...using the CERN algorithm.... 
cwerf            - Marco representation
complexErrFunc   - Serge complex representation
*/

double VoigtA(const double xi, const double x0, const double a, const double wD, const double w0, const double d0)
{
	double rx = (M_SqrtLn2 * (xi - (x0 + d0)) / wD);
	double ry = ((w0 / wD) * M_SqrtLn2);
	double vr, vi, vr0;
	cwerf(0.0, ry, &vr0, &vi);
	cwerf(rx, ry, &vr, &vi);
	return a * vr / vr0;
}

double VoigtS(const double xi, const double x0, const double s, const double wD, const double w0, const double d0)
{
	dcplx z((M_SqrtLn2 * (xi - (x0 + d0)) / wD), ((w0 / wD) * M_SqrtLn2));

	double rx = (M_SqrtLn2 * (xi - (x0 + d0)) / wD);
	double ry = ((w0 / wD) * M_SqrtLn2);
	double vr, vi;
	cwerf(rx, ry, &vr, &vi);
	return s * vr * (1 / wD) * M_SqrtLn2oPi;
}

/*
VoigtH peak

Vs(x; wD,w0) = SQRT(ln2/PI) * Re[w(z)]/wD
z = SQRT(ln2) * (x/wD) , SQRT(ln2) * (w0/wD)
...the w is calculated with the Humlicek algorithm
*/

double VoigtHA(const double xi, const double x0, const double a, const double wD, const double w0, const double d0)
{
	dcplx z((M_SqrtLn2 * (xi - (x0 + d0)) / wD), ((w0 / wD) * M_SqrtLn2));
	dcplx z0(0.0, ((w0 / wD) * M_SqrtLn2));
	//return a * CPF(z).real() / CPF(z0).real();
	return a * (Cerf::faddeeva_fast(z).real()) / (Cerf::faddeeva(z0).real());
}

double VoigtHS(const double xi, const double x0, const double s, const double wD, const double w0, const double d0)
{
	dcplx z((M_SqrtLn2 * (xi - (x0 + d0)) / wD), ((w0 / wD) * M_SqrtLn2));
	//return s * (CPF(z).real()) * (1 / wD) * M_SqrtLn2oPi;
	  return s * (Cerf::faddeeva_fast(z).real()) * (1 / wD) * M_SqrtLn2oPi;
}


/*
the HTP function...
C0 = w0 + i d0
C2 = w2 + i d2

...constants according to the papier of Forthomme...


*/
double HTPA(double &xi, double *Par)
{

	//some constants for complex
	static const dcplx iCplx((double)0.0, (double)1.0);
	static const dcplx ZEROCplx((double)0.0, (double)0.0);

	double realIN(0.0);
	double realOUT(0.0);
	double imagIN(0.0);
	double imagOUT(0.0);

	//variables to determine a good region....
	long region(0);
	double SZp(0.0);
	double SZm(0.0);
	double DSZ(0.0);
	double SZmax(0.0);
	double SZmin(0.0);

	const double nuVC = abs(*(Par + 10));
	const double eta = abs(*(Par + 9));

	const double va0 = P_CoSqrtLn2 *  ((*(Par + 4)) / (*(Par + 1))); // va0 = (c/sqrt(ln2)) * wD/nu0
	const double x0c(P_C / (va0 * (*(Par + 1))));

	const dcplx C0{ *(Par + 5), 1.0 * (*(Par + 6)) };
	const dcplx C2{ *(Par + 7), 1.0 * (*(Par + 8)) };
	
	const dcplx C0Tilde(((double)1.0 - eta) * (C0 - (double)1.5 * C2) + nuVC);
	const dcplx C2Tilde(C2 * ((double)1.0 - eta));
	
	const dcplx ix  (0.0, (*(Par + 1)) - xi);
	
	dcplx invC2Tilde(0.0, 0.0);
	dcplx X(0.0, 0.0);
	dcplx Yr(0.0, 0.0);
	dcplx Y(0.0, 0.0);
	
	dcplx Zp(0.0, 0.0);
	dcplx Zm(0.0, 0.0);
	
	dcplx wZp(0.0, 0.0);
	dcplx wZm(0.0, 0.0);
	
	dcplx Aterm(0.0, 0.0);
	dcplx Bterm(0.0, 0.0);

	dcplx qpcSDNGP(0.0, 0.0);

	// selecting regions...
	if (C2Tilde != ZEROCplx) 
	{
		invC2Tilde = 1.0 / C2Tilde;
		X = ( (ix + C0Tilde) * invC2Tilde );
		Yr = (*(Par + 1)) * va0 * invC2Tilde * P_1o2C;
		Y = Yr*Yr;
	}
	else
	{
	// when C2t=0.0
		region = 3;
	}
	// when abs(X) is much larger than abs(Y)
	if (abs(Y)<(double)1.0e-15*abs(X))
		region = 1;
	// when abs(Y) is much larger than abs(X)
	else if (abs(X)<(double)3.0e-8*abs(Y))
		region = 2;

	switch (region)
	{
	
	case 0: {
		invC2Tilde = 1.0 / C2Tilde;
		// calculating Z1 and Z2
		Zm = sqrt(X + Y) - Yr;
		Zp = Zm + 2.0*Yr;
		// calculating the real and imaginary parts of Z1 and Z2
		//..xz1 = -z1CPlx.imag();
		//..yz1 = z1CPlx.real();
		//..xz2 = -z2CPlx.imag();
		//..yz2 = z2CPlx.real();
		// check if Z1 and Z2 are close to each other ... well we ca use the abs function here nope???
		SZp = abs(Zp);//..SZ2 = fcnNorme(xz2, yz2);
		SZm = abs(Zm);//..SZ1 = fcnNorme(xz1, yz1);
		DSZ = abs(SZm - SZp);
		SZmax = (((SZm) >= (SZp)) ? (SZm) : (SZp));  //  max(SZm, SZp);
		SZmin = (((SZm) <= (SZp)) ? (SZm) : (SZp));  //  min...
		// when Z1 and Z2 are close to each other, ensure that they are in
		// the same interval of CPF
		if (DSZ <= (double)1.0 && SZmax > (double)8.0 &&  SZmin <= (double)8.0)
		{
			wZp = CPF3(iCplx * Zp);
			wZm = CPF3(iCplx * Zm);
		}
		else
		{
			wZp = CPF(iCplx * Zp);
			wZm = CPF(iCplx * Zm);
		}
		// calculating the A term of the profile
		//..Aterm = *(ptrParamP)*SQRTPI*(w1CPlx - w2CPlx);
		//..Bterm = (-(double)1.0 + SQRTPI / ((double)2.0*yCplx)*((double)1.0 - fcnCarre(z1CPlx))*w1CPlx
		//..	                  - SQRTPI / ((double)2.0*yCplx)*((double)1.0 - fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
		
		Aterm = M_SqrtPI * x0c * (wZm - wZp);
		Bterm =  invC2Tilde * va0 * va0 * (  -(double)1.0 
			+ M_SqrtPI / ((double)2.0*Yr) * ((double)1.0 - Zm*Zm) * wZm 
			- M_SqrtPI / ((double)2.0*Yr) * ((double)1.0 - Zp*Zp) * wZp );
	}
			break;
	///////////////////////////////////////////////////////////////////
	case 1: {
		invC2Tilde = 1.0 / C2Tilde;
		// we need SQRT(X) and SQRT(X+Y) only
		Zm = sqrt(X);     
		Zp = sqrt(X + Y); 
		// the error func. are always acompaigned by ther arguments
		wZm = Zm * CPF(iCplx * Zm);
		wZp = Zp * CPF(iCplx * Zp);

		// when abs(X) is much larger than abs(Y)
		if (abs(sqrt(X)) <= (double)4000.0)
		{
			Aterm = ((double)2.0*M_SqrtPI*invC2Tilde)*(M_1oSqrtPI - wZm);
			Bterm = invC2Tilde * va0 * va0 * (-(double)1.0 + (double)2.0*M_SqrtPI*((double)1.0 - X - (double)2.0*Y)*(M_1oSqrtPI - wZm)
				+ (double)2.0*M_SqrtPI*wZp);
		}
		else
			// when abs(X) is much larger than 1
		{
			Aterm = invC2Tilde*((double)1.0/X - (double)1.5/(X*X));
			Bterm = invC2Tilde*(-(double)1.0 + ((double)1.0 - X - (double)2.0*Y)*((1.0/X) - (double)1.5*(1.0/(X*X)))
				+ (double)2.0*M_SqrtPI*wZp);
		}
	}
			break;
	///////////////////////////////////////////////////////////////////
	case 2: {
		invC2Tilde = 1.0 / C2Tilde;
		// when abs(Y) is much larger than abs(X)
		//... dcplx c0t(*(ptrParamP + 1), *(ptrParamP + 2)); already there
		Zm = (ix + C0Tilde) * x0c;//..z1CPlx = *(ptrParamP)*(iCplx*x + c0t);
		Zp = sqrt(X + Y) + Yr;//..z2CPlx = sqrt(xCplx + y2Cplx) + yCplx;
		//..xz1 = -z1CPlx.imag();
		//..yz1 = z1CPlx.real();
		//..xz2 = -z2CPlx.imag();
		//..yz2 = z2CPlx.real();    ...perhaps similar case for CPF calc. should apply here?
		wZm = CPF(iCplx * Zm);          //w1CPlx = CPF(dcplx(xz1, yz1));
		wZp = CPF(iCplx * Zp);          //w2CPlx = CPF(dcplx(xz2, yz2));
		
		Aterm = M_SqrtPI * x0c * (wZm - wZp); 
		//Aterm = *(ptrParamP)*SQRTPI*(w1CPlx - w2CPlx);
		Bterm = invC2Tilde * va0 * va0 * (-(double)1.0
			+ M_SqrtPI / ((double)2.0*Yr) * ((double)1.0 - Zm*Zm) * wZm
			- M_SqrtPI / ((double)2.0*Yr) * ((double)1.0 - Zp*Zp) * wZp);
		//Bterm = (-(double)1.0 + SQRTPI / ((double)2.0*yCplx)*((double)1.0 - fcnCarre(z1CPlx))*w1CPlx
		//	                    - SQRTPI / ((double)2.0*yCplx)*((double)1.0 - fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
	}
			break;
	///////////////////////////////////////////////////////////////////
	case 3: {
		// when C2t=0
		Zm = ((ix + C0Tilde) * x0c);                            //eq. 8-1
		wZm = Cerf::faddeeva(iCplx * Zm);                              
		Aterm = M_SqrtPI * x0c * wZm;
		if (abs(sqrt(Zm))<(double)4000.0)
		{
			Bterm = M_SqrtPI * x0c * va0 * va0 *  (((1.0 - (Zm*Zm)) * wZm) + Zm * M_1oSqrtPI);  // Bterm = *(ptrParamP)*SQRTPI*((((double)1.0-fcnCarre(z1CPlx))*w1CPlx) + *(ptrParamP)*z1CPlx);
		}
		else
		{
			Bterm = M_SqrtPI * x0c * va0 * va0 *  (wZm + ((double)0.5 / Zm) - ((double)0.75 *  M_1oSqrtPI * (Zm * Zm *Zm))); //Bterm = *(ptrParamP)*((SQRTPI*w1CPlx) + (double)0.5 / z1CPlx - (double)0.75*Power<3>::of(z1CPlx));
		}
	}
			break;
	default:break;
	};

	qpcSDNGP = M_1_PI * (Aterm/(1.0 - (nuVC-eta*(C0-1.5*C2))*Aterm + (eta*C2/(va0*va0))*Bterm)); //( ((double)1.0 - nuVC - eta*(C0 - (double)1.5*C2))*Aterm + (eta * C2 /(va0 * va0))*Bterm)    );
	return ((*(Par + 2))*qpcSDNGP.real());//((*(Par + 2))*qpcSDNGP.real());
}
















// ************************************************************************ //
/*         Fonction d'Erreur dcplx w(z): Fonction FADDEEVA                  */
// ************************************************************************ //

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

/*********/
/* cwerf */
/*********/

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

dcplx complexErrFunc(const dcplx &z)
{
	// ======================================================================================== //
	/* Return CERNlib complex error function                                                    */
	/*                                                                                          */
	/* This code is translated from the fortran version in the CERN mathlib.                    */
	/* (see http://aliceinfo.cern.ch/alicvs/viewvc/MINICERN/mathlib/gen/c/cwerf64.F?view=co)    */
	// ======================================================================================== //

	static concurrent_vector<dcplx> r(38);
	// { HF, C1, C2,    C3, C4,                                P}
	static double lookup_table[] = { 0.5,7.4,8.3,0.3125,1.6,46768052394588893.382517914646921 };
	dcplx zh, s, t, v;
	static const dcplx zero(0.0, 0.0);
	double x(z.real()), y(z.imag()), xAbs(fabs(x)), yAbs(fabs(y));
	int N;

	if ((yAbs < lookup_table[1]) && (xAbs <  lookup_table[2]))
	{
		zh = dcplx(yAbs + lookup_table[4], xAbs);
		r[37] = zero;
		N = 36;
		while (N > 0)
		{
			t = zh + conj(r[N + 1])*(double)N;
			r[N--] = (t*lookup_table[0]) / norm(t);
		}
		double xl = lookup_table[5];
		s = zero;
		N = 33;
		while (N > 0)
		{
			xl = lookup_table[3] * xl;
			s = r[N--] * (s + xl);
		}
		v = s * M_2SqrtPI;
	}
	else
	{
		zh = dcplx(yAbs, xAbs);
		r[1] = zero;
		N = 9;
		while (N > 0)
		{
			t = zh + conj(r[1])*(double)N;
			r[1] = (t*lookup_table[0]) / norm(t);
			N--;
		}
		v = r[1] * M_2SqrtPI;
	}
	if (yAbs == (double)0.0)
		v = dcplx(exp(-(xAbs*xAbs)), v.imag());
	if (y < (double)0.0)
	{
		dcplx tmp(yAbs*yAbs - xAbs*xAbs, -(double)2.0*xAbs*yAbs);
		v = (double)2.0*exp(tmp) - v;
		if (x >(double)0.0)
			v = conj(v);
	}
	else
	{
		if (x < (double)0.0)
			v = conj(v);
	}
	return v;
}


/*Marco representation.......
from CERN mathlib cwerf64.F
.........*/

void cwerf(double x, double y, double *vr, double *vi)									//*
{																						//*
																						//*
	int n;																				//*
	double xa, ya;																			//*
	double zhr, zhi;																		//*
	double rr[38], ri[38];																	//*
	double tr, ti;																			//*
	double r, xl;																			//*
	double sr, si;																			//*
																							//*
	xa = fabs(x);																			//*
	ya = fabs(y);																			//*
	if (ya<7.4 && xa<8.3)																	//*
	{																					//*	
		zhr = ya + 1.6;																		//*
		zhi = xa;																			//*
		rr[37] = 0;																			//*
		ri[37] = 0;																			//*
		for (n = 36; n>0; n--)																//*
		{																					//*
			tr = zhr + n*rr[n + 1];																	//*
			ti = zhi - n*ri[n + 1];																	//*
			xl = 0.5 / (tr*tr + ti*ti);																//*
			rr[n] = tr*xl;																		//*
			ri[n] = ti*xl;																		//*
		}																					//*
		xl = 46768052394589056.;															//*
		sr = 0;																				//*
		si = 0;																				//*
		for (n = 33; n>0; n--)																//*
		{																					//*
			xl *= .3125;																		//*
			tr = rr[n] * (sr + xl) - ri[n] * si;														//*
			si = rr[n] * si + ri[n] * (sr + xl);														//*
			sr = tr;																			//*
		}																					//*
		*vr = sr*CVOIGT;																	//*
		*vi = si*CVOIGT;																	//*
	}																					//*
	else																					//*
	{																					//*				
		zhr = ya;																			//*
		zhi = xa;																			//*
		rr[0] = 0.;																			//*
		ri[0] = 0.;																			//*
		for (n = 9; n>0; n--)																	//*		
		{																					//*
			tr = zhr + n*rr[0];																	//*
			ti = zhi - n*ri[0];																	//*
			xl = 0.5 / (tr*tr + ti*ti);																//*
			rr[0] = tr*xl;																		//*
			ri[0] = ti*xl;																		//*
		}																					//*
		*vr = rr[0] * CVOIGT;																	//*
		*vi = ri[0] * CVOIGT;																	//*
	}																					//*
	if (ya == 0.)																			//*
		*vr = exp(-xa*xa);																	//*
	if (y<0.)																				//*	
	{																					//*
		tr = xa*xa - ya*ya;																	//*
		ti = 2.*xa*ya;																		//*
		r = sqrt(tr*tr + ti*ti);																//*
		xl = exp(-r);																		//*
		sr = 2.*xl*tr / r;																	//*
		si = -2.*xl*ti / r;																	//*
		*vr = 2.*sr - *vr;																	//*	
		*vi = 2.*si - *vi;																	//*
		if (x>0.)																			//*
			*vi = -*vi;																			//*		
	}																					//*
	else																					//*
	{																					//*
		if (x<0)																			//*
			*vi = -*vi;																			//*
	}																					//*	
	return;																				//*
}



// ============================================================================================================================================================================================================================ //
//																									"CPF": Complex Probability Function
// ============================================================================================================================================================================================================================ //
//         .................................................
//         .       Subroutine to Compute the Complex       .
//         .        Probability Function W(z=X+iY)         .
//         .    W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0       .
//         .    Which Appears when Convoluting a Complex   .
//         .    Lorentzian Profile by a Gaussian Shape     .
//         .................................................
//
//
// This Routine was Taken from the Paper by J. Humlicek, which
// is Available in Page 309 of Volume 21 of the 1979 Issue of
// the Journal of Quantitative Spectroscopy and Radiative Transfer
// Please Refer to this Paper for More Information
//
// Accessed Files:  None
// --------------
//
// Called Routines: None
// ---------------
//
// ============================================================================================================================================================================================================================ //

dcplx CPF(const dcplx &z)
{

	dcplx  result(0.0, 0.0);
	double x(z.real());
	double y(z.imag());
	const double T[] = {
		0.314240376,
		0.947788391,
		1.59768264,
		2.27950708,
		3.02063703,
		3.8897249
	};
	const double U[] = {
		1.01172805,
		-0.75197147,
		1.2557727e-2,
		1.00220082e-2,
		-2.42068135e-4,
		5.00848061e-7,
	};
	const double S[] = {
		1.393237,
		0.231152406,
		-0.155351466,
		6.21836624e-3,
		9.19082986e-5,
		-6.27525958e-7
	};
	long region(2);

	// new Region 3
	if (abs(z)>(double)8.0)
		region = 3;
	else if (y>(double)0.85 || fabs(x)<((double)18.1*y + (double)1.65))
		region = 1;
	switch (region)
	{
	case 1: {
		const double Y1 = y + (double)1.5;
		const double Y2 = Y1*Y1;
		double wr(0.0);
		double wi(0.0);
		double R(0.0);
		double R2(0.0);
		double D(0.0);
		double D1(0.0);
		double D2(0.0);
		double D3(0.0);
		double D4(0.0);

		for (long i = 0; i<6; ++i)
		{
			R = x - T[i];
			R2 = R*R;
			D = (double)1.0 / (R2 + Y2);
			D1 = Y1*D;
			D2 = R*D;
			R = x + T[i];
			R2 = R*R;
			D = (double)1.0 / (R2 + Y2);
			D3 = Y1*D;
			D4 = R*D;
			wr += (U[i] * (D1 + D3)) - (S[i] * (D2 - D4));
			wi += (U[i] * (D2 + D4)) + (S[i] * (D1 - D3));
		}
		result = dcplx(wr, wi);
	}
			break;
	case 2: {
		const double Y1 = y + (double)1.5;
		const double Y2 = Y1*Y1;
		const double Y3 = y + (double)3.0;
		double wr(0.0);
		double wi(0.0);
		double R(0.0);
		double R2(0.0);
		double D(0.0);
		double D1(0.0);
		double D2(0.0);
		double D3(0.0);
		double D4(0.0);

		if (abs(x)<(double)12.0)
			wr = exp(-x*x);
		for (long i = 0; i<6; ++i)
		{
			R = x - T[i];
			R2 = R*R;
			D = (double)1.0 / (R2 + Y2);
			D1 = Y1*D;
			D2 = R*D;
			wr += y*(U[i] * (R*D2 - (double)1.5*D1) + S[i] * Y3*D2) / (R2 + (double)2.25);
			R = x + T[i];
			R2 = R*R;
			D = (double)1.0 / (R2 + Y2);
			D3 = Y1*D;
			D4 = R*D;
			wr += y*(U[i] * (R*D4 - (double)1.5*D3) - S[i] * Y3*D4) / (R2 + (double)2.25);
			wi += U[i] * (D2 + D4) + S[i] * (D1 - D3);
		}
		result = dcplx(wr, wi);
	}
			break;
	case 3: {
		const dcplx zone(1.0, 0.0);
		const dcplx zi(0.0, 1.0);
		const double tt[] = {
			0.5,
			1.5,
			2.5,
			3.5,
			4.5,
			5.5,
			6.5,
			7.5,
			8.5,
			9.5,
			10.5,
			11.5,
			12.5,
			13.5,
			14.5
		};
		dcplx zm1(0.0, 0.0);
		dcplx zm2(0.0, 0.0);
		dcplx zterm(0.0, 0.0);
		dcplx zsum(0.0, 0.0);

		zm1 = zone / z;
		zm2 = zm1*zm1;
		zsum = zone;
		zterm = zone;
		for (long i = 0; i<15; ++i)
		{
			zterm *= zm2*tt[i];
			zsum += zterm;
		}
		zsum *= zi*zm1* M_1oSqrtPI;
		result = zsum;
	}
			break;
	}

	return result;
}

// ============================================================================================================================================================================================================================ //
//																									"CPF3": Complex Probability Function
// ============================================================================================================================================================================================================================ //
//         .................................................
//         .       Subroutine to Compute the Complex       .
//         .        Probability Function W(z=X+iY)         .
//         .    W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0       .
//         .    Which Appears when Convoluting a Complex   .
//         .    Lorentzian Profile by a Gaussian Shape     .
//         .................................................
//
//
// This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8
//
// Accessed Files:  None
// --------------
//
// Called Routines: None
// ---------------
//
// ============================================================================================================================================================================================================================ //

dcplx CPF3(const dcplx &z)
{
	dcplx zone(1.0, 0.0);
	dcplx zi(0.0, 1.0);
	dcplx zm1(0.0, 0.0);
	dcplx zm2(0.0, 0.0);
	dcplx zterm(0.0, 0.0);
	dcplx zsum(0.0, 0.0);
	double pipwoeronehalf(0.564189583547756);
	double tt[] = {
		0.5,
		1.5,
		2.5,
		3.5,
		4.5,
		5.5,
		9.5,
		10.5,
		11.5,
		12.5,
		13.5,
		14.5
	};

	zm1 = zone / z;
	zm2 = zm1*zm1;
	zsum = zone;
	zterm = zone;
	for (int i = 0; i<15; ++i)
	{
		zterm *= zm2*tt[i];
		zsum += zterm;
	}
	zsum *= zi*zm1*pipwoeronehalf;

	return zsum;
}





