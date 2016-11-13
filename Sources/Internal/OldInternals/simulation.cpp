//http://bytes.com/topic/c/answers/849132-std-vector-c-array
//memcpy vector c++

#include "simulation.h"


// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

namespace ROOT {
	namespace Minuit2 {

// ============================================================================================================================================================================================================================ //
//								                                         Typedef                                                                                                                                                //
// ============================================================================================================================================================================================================================ //

// Pour les différents types de profil
typedef  long (simulation::*ptrDoProfile)(long &, double *, std::vector<double>::const_iterator);
// Pour le calcul des surfaces des différents profils
typedef  double (simulation::*ptrDoAreaProfile)(std::vector<double>::const_iterator);

// ============================================================================================================================================================================================================================ //
//																							Constructeurs											                                                                            //
// ============================================================================================================================================================================================================================ //

	// constructeur standard
	simulation::simulation(const std::vector<double> &x,const std::vector<double> &par):fX(x),fPar(par) {}

	// constructeur pour initialiser les paramètres uniquement
	simulation::simulation(const std::vector<double> &par):fPar(par)
	{
		 fX = std::vector<double>();
	}

	// constructeur par défaut
	simulation::simulation()
	{
		fX   = std::vector<double>();
		fPar = std::vector<double>();
	}


// ============================================================================================================================================================================================================================ //
//																				 Destructeur																																	//
// ============================================================================================================================================================================================================================ //

	simulation::~simulation() {}

// ============================================================================================================================================================================================================================ //



// ============================================================================================================================================================================================================================ //
//								                                                   Mutateurs                                                                                                                                    //
// ============================================================================================================================================================================================================================ //

void simulation::setParam(const long &index, const double &newValue)
{
	fPar[index]=newValue;
}

void simulation::setParamProfile(const std::vector<double> &x,const std::vector<double> &par)
{
	if (fPar.size() != par.size())
		fPar.resize(par.size());
	if (fX.size() != x.size())
		fX.resize(x.size());
	fPar = par;
	fX = x;
}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

									/******************************************************************/
									/*  Fonction d'appel des profils pour le calcul d'une simulation  */
									/******************************************************************/

// ============================================================================================================================================================================================================================ //


	// Par construction le vecteur par contient le nombre de courbes,
	// le nombre de coefficients du polynome simulant la ligne de base
	// et si méthode liée ou pas
	// par[0] = #lines
	// par[1] = #coefficients du polynome
	// par[2] = Amplitude (=0) or area (=1)
	// par[3] = method
	// standard method = 0
	// fit parameter-by-parameter = 1
	// Linked parameters
	// Pressure

void simulation::profileSimulation(double *y)
{
	long typeSimulation=(long)fPar[3];
	switch (typeSimulation)
	{
		case 0: // STANDARD
				{
					stdSimulationProfile(y);
				}
				break;
				// PARAMATER BY PARAMETER
		case 1: break;
				// LINKED
		case 2: break;
				// PRESSURE
		case 3: break;
		default: break;
	}
}


/***********************************/
/* simulation dans le cas standard */
/***********************************/

void simulation::stdSimulationProfile(double *y)
{

	long nbLines = (long)fPar[0];
	long nbCoeffsPolynome = (long)fPar[1];
	long area = (long)fPar[2];
	long line(0);
	long typeSimul(2);
	long res(0);
	std::vector<double>::const_iterator ptrParam = fPar.begin()+ARRAYS_SIZE;
	std::vector<double>::const_iterator ptrParamEnd = ptrParam + nbLines*NB_PARAM_STD;
    static ptrDoProfile tabDoProfile[] = {
									       &simulation::doDoppler,
										   &simulation::doLorentz,
		                                   &simulation::doVoigt,
										   &simulation::doHumlicek,
										   &simulation::doNelkinGhatak,
										   &simulation::doNelkinGhatakHumlicek,
										   &simulation::doGalatry,
										   &simulation::doIncompleteGamma,
										   &simulation::doFano,
										   &simulation::doQuadraticSpeedDependentVoigt,										   
										   &simulation::doQuadraticSpeedDependentNelkinGhatak,
										   &simulation::doQuadraticSpeedDependentVoigtBoone,
										   &simulation::doQuadraticSpeedDependentNelkinGhatakBoone,
										   &simulation::doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone,
										   &simulation::doQuadraticSpeedDependentHumlicekBoone,
										   &simulation::doQuadraticSpeedDependentNelkinGhatakHumlicekBoone,
										   &simulation::doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone,
										   &simulation::doApproximateQuadraticSpeedDependentGalatry
		                                };

	if (nbCoeffsPolynome)
		FCN_Baseline(y,fPar);

	if (nbLines != 1)
	{
		long beginNbProfile(0);
		long endNbProfile(nbLines - 1);
		long nbProfile(0);
		long indexParam(0);
		long indexParamPlus1(0);
		bool nextIteration(true);

		while (nextIteration && nbProfile<endNbProfile)
		{
			indexParam = ARRAYS_SIZE + nbProfile*NB_PARAM_STD;
			indexParamPlus1 = ARRAYS_SIZE + (nbProfile+1)*NB_PARAM_STD;
			if (fPar[indexParam]!=fPar[indexParamPlus1])
				nextIteration = false;
			++nbProfile;
		}
		typeSimul = 0;
		if (nextIteration == false)
			typeSimul = 1;
	}
	switch(typeSimul)
	{
		case 0:{
					std::vector<double>::const_iterator ptrCourant;
					line = (long)(*ptrParam);
					for (long i = 0;i<nbLines;++i)
					{
						ptrCourant = i*NB_PARAM_STD + ptrParam;
						res = (*this.*tabDoProfile[line])(area,y,ptrCourant);
					}
			   }
			   break;
		case 1:{
					std::vector<double>::const_iterator ptrCourant;
					for (long i = 0;i<nbLines;++i)
					{
						ptrCourant = i*NB_PARAM_STD + ptrParam;
						line = (long)(*ptrCourant);
						res = (*this.*tabDoProfile[line])(area,y,ptrCourant);
					}
			   }
			   break;
		case 2:{
					line = (long)(*ptrParam);
					res = (*this.*tabDoProfile[line])(area,y,ptrParam);
			   }
			   break;
	}
}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

															/********************************************************/
															/* Fonction d'appel pour le calcul des aires du spectre */
															/* La normalisation est: 1.0/f(x=0)                     */
															/********************************************************/

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

void simulation::areaProfile(double *profileArea)
{
	long nbLines = (long)fPar[0];
	long typeSimul(2);
	long line(0);
	std::vector<double>::const_iterator ptrParam = fPar.begin()+ARRAYS_SIZE;
	std::vector<double>::const_iterator ptrParamEnd = ptrParam + nbLines*NB_PARAM_STD;
	ptrDoAreaProfile tabDoAreaProfile[] = {
										   &simulation::doAreaDoppler,
										   &simulation::doAreaLorentz,
		                                   &simulation::doAreaVoigt,
										   &simulation::doAreaHumlicek,
										   &simulation::doAreaNelkinGhatak,
										   &simulation::doAreaNelkinGhatakHumlicek,
										   &simulation::doAreaGalatry,
										   &simulation::doAreaIncompleteGamma,
										   &simulation::doAreaFano,
										   &simulation::doAreaQuadraticSpeedDependentVoigt,
										   &simulation::doAreaQuadraticSpeedDependentNelkinGhatak,
										   &simulation::doAreaQuadraticSpeedDependentVoigtBoone,
										   &simulation::doAreaQuadraticSpeedDependentNelkinGhatakBoone,
										   &simulation::doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone,
										   &simulation::doAreaQuadraticSpeedDependentHumlicekBoone,
										   &simulation::doAreaQuadraticSpeedDependentNelkinGhatakHumlicekBoone,
										   &simulation::doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone,
										   &simulation::doAreaApproximateQuadraticSpeedDependentGalatry
		                                  };

	if (nbLines != 1)
	{
		long beginNbProfile(0);
		long endNbProfile(nbLines - 1);
		long nbProfile(0);
		long indexParam(0);
		long indexParamPlus1(0);
		bool nextIteration(true);

		while (nextIteration && nbProfile<endNbProfile)
		{
			indexParam = ARRAYS_SIZE + nbProfile*NB_PARAM_STD;
			indexParamPlus1 = ARRAYS_SIZE + (nbProfile+1)*NB_PARAM_STD;
			if (fPar[indexParam]!=fPar[indexParamPlus1])
				nextIteration = false;
			++nbProfile;
		}
		typeSimul = 0;
		if (nextIteration == false)
			typeSimul = 1;
	}
	switch(typeSimul)
	{
		case 0:{
					std::vector<double>::const_iterator ptrCourant;
					line = (long)(*ptrParam);
					for (long i = 0;i<nbLines;++i)
					{
						ptrCourant = i*NB_PARAM_STD + ptrParam;
						profileArea[i] = (*this.*tabDoAreaProfile[line])(ptrCourant);
					}
			   }
			   break;
		case 1:{
					std::vector<double>::const_iterator ptrCourant;
					for (long i = 0;i<nbLines;++i)
					{
						ptrCourant = i*NB_PARAM_STD + ptrParam;
						line = (long)(*ptrCourant);
						profileArea[i] = (*this.*tabDoAreaProfile[line])(ptrCourant);
					}
			   }
			   break;
		case 2:{
					line = (long)(*ptrParam);
					*profileArea = (*this.*tabDoAreaProfile[line])(ptrParam);
			   }
	}


}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																							/***********************************/
																							/* SPECTRAL LINE SHAPES SIMULATION */
																							/***********************************/

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //


// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																							/*****************************************/
																							/* BASIC MODELS FOR SPECTRAL LINE SHAPES */
																							/*****************************************/

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //


// ============================================================================================================================================================================================================================ //

																								//****************************//
																								// DOPPLER PROFILE SIMULATION //
																								//****************************//

// ============================================================================================================================================================================================================================ //

long simulation::doDoppler(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	std::vector<double> paramProfile(1,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double center = *(ptrParam+Wavelength);
	const double invGWidth = fcnInverse(gHWHM);
	const long typeSimul = (!(*(ptrParam+parameter_0)) ? 0:1);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;

	const double amp = (!area) ? (amplitude/dopplerProfile(ptrParamProfile,lineCenter)):amplitude;

	switch (typeSimul)
	{
		case 0:{
					for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;
						(*ptrY++)+=amp*dopplerProfile(ptrParamProfile,lineCenter);
					}			   
			   }
			   break;
		case 1:{
					double sigmaDoppler = fabs(*(ptrParam+parameter_0));
					for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx,++ptrY)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;
						if (fabs(lineCenter)<=sigmaDoppler)
							(*ptrY)+=amp*dopplerProfile(ptrParamProfile,lineCenter);
						else if (lineCenter>=sigmaDoppler)
							break;
					}			   
			   }
			   break;
		default: break;
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//*****************//
																						// DOPPLER PROFILE //
																						//*****************//

// ============================================================================================================================================================================================================================ //

double simulation::dopplerProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	double lineCenter2 = fcnCarre(x);

	return ((*ptrParamP)*exp(-lineCenter2));
}

// ============================================================================================================================================================================================================================ //

																						//**********************//
																						// AREA DOPPLER PROFILE //
																						//**********************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaDoppler(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return dblERREUR;
	if(!amplitude)
		return dblOK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblOK;

	const double norm = SQRTPI*gHWHM;

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																						//****************************//
																						// LORENTZ PROFILE SIMULATION //
																						//****************************//

// ============================================================================================================================================================================================================================ //

long simulation::doLorentz(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if(!lHWHM)
		return OK;

	std::vector<double> paramProfile(2,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double center = *(ptrParam+Wavelength);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_PI;
	paramProfile[1] = lHWHM;

	const double amp = (!area) ? (amplitude/lorentzProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
			lineCenter = ((*ptrFx)-center);
			(*ptrY++)+=amp*lorentzProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//*****************//
																						// LORENTZ PROFILE //
																						//*****************//

// ============================================================================================================================================================================================================================ //

double simulation::lorentzProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	double lineCenter2 = fcnCarre(x);
	double gamma2 = fcnCarre(*(ptrParamP+1));
	double intensity = *(ptrParamP)*(*(ptrParamP+1))*fcnInverse(gamma2+lineCenter2);

	return intensity;
}

// ============================================================================================================================================================================================================================ //

																						//**********************//
																						// AREA LORENTZ PROFILE //
																						//**********************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaLorentz(std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return dblERREUR;
	if(!amplitude)
		return dblOK;

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if(!lHWHM)
		return dblOK;

	double norm = PI*lHWHM;

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //
//
// For simplicity and convenience, dimensionless parameters are introduced following Herbert (a)
// and Varghese and Hanson (b):
//		xTilde = (W-Wfi)/DWd (ie Wfi=W0)
//		x = (W-Wfi-delta)/DWd = xTilde-s
//		with
//		s = Delta/DWd
//		y = Gamma/DWd
//
// where gamma*P=Gamma is the collisional HWHM and d*P=Delta is the collisional shift at pressure P.
//
// For this simulation, the collisional shifting Delta is set to zero
//
// DWd (in cm-1), which is the 1/e Doppler half width, different from the more commonly used
// half width at half maximum denoted gammaD, expressed as:
//
//	gammaD = sqrt(ln2)*DWd/2*PI = Wfi*sqrt(((ln2)*kb*T)/(2*PI*PI*m)) = 3.58*10-7*Sigmafi*sqrt(T/M)
//
// where:
// 	 * M is the molar mass in g/mol;
//	 * the wavenumber Sigmafi in cm-1;
//	 * Wfi is the angular frequency (in rad/s) for the considered f<-i
//	  optical transition at zero pressure (ie in the absence of any collisional effects).
//
// ============================================================================================================================================================================================================================ //
//.
// (a) Herbert F. Spectrum line profiles: a generalized Voigt function including collisional narrowing.
//    JQSRT 14, 943 (1974)
// (b) Varghese PL, Hanson RK. Collisional narrowing effects on spectral line shapes measured at high resolution.
//    Appl Opt 23, 2376 (1984)
//
// ============================================================================================================================================================================================================================ //

																						//**************************//
																						// VOIGT PROFILE SIMULATION //
																						//**************************//

// ============================================================================================================================================================================================================================ //
// In this case, the  complex  probability function w(z) is computed thanks to the CERNlib routine
// ============================================================================================================================================================================================================================ //
// Usage:
//	In any arithmetic expression, CWERF(Z) has the value w(Z) ,
//  where CWERF is of type COMPLEX, WWERF is of type COMPLEX*16, and Z has the same type as the function name.
// Method:
//  The method is described in Ref. 2.
// Accuracy:
//  CWERF (except on CDC and Cray computers) has full single-precision accuracy. 
//  For most values of the argument Z, WWERF (and CWERF on CDC and Cray computers) has an accuracy of approximately two significant digits 
//  less than the machine precision.
// Notes:
//  This subprogram is a modified version of the algorithm presented in Ref. 1.
// ============================================================================================================================================================================================================================ //
// References:
//1. W. Gautschi, Algorithm 363, Complex Error Function, Collected Algorithms from CACM (1969).
//2. W. Gautschi, Efficient Computation of the Complex Error Function, SIAM J. Numer. Anal. 7 (1970) 187--198.
//3. K.S. Kölbig, Certification of Algorithm 363 Complex Error Function, Comm. ACM 15 (1972) 465--466.
// ============================================================================================================================================================================================================================ //

long simulation::doVoigt(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return doLorentz(area,y,ptrParam);

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if (!lHWHM)
		return doDoppler(area,y,ptrParam);

	std::vector<double> paramProfile(2,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
    const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;
	const double center = *(ptrParam+Wavelength) + delta;
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = fabs(lHWHM*invGWidth);

	const double amp = (!area) ? (amplitude/voigtProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*voigtProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//***************//
																						// VOIGT PROFILE //
																						//***************//

// ============================================================================================================================================================================================================================ //

double simulation::voigtProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	dcplx z(x,(*(ptrParamP+1)));
	dcplx cwerf = complexErrFunc(z) ;
	double intensity = (*ptrParamP)*cwerf.real();

    return intensity;
}

// ============================================================================================================================================================================================================================ //

																						//********************//
																						// AREA VOIGT PROFILE //
																						//********************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaVoigt(std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return dblERREUR;
	if(!amplitude)
		return dblOK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return doAreaLorentz(ptrParam);

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if (!lHWHM)
		return doAreaDoppler(ptrParam);

	const double invGWidth = abs(fcnInverse(gHWHM));
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(2,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = fabs((*(ptrParam+gammaHWHM))*invGWidth);
	const double norm = fcnInverse(voigtProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}
// ============================================================================================================================================================================================================================ //

																						//***************************************************//
																						// VOIGT PROFILE SIMULATION (HUMLICEK APPROXIMATION) //
																						//***************************************************//

// ============================================================================================================================================================================================================================ //
// Based on a slightly improved version of the CPF subroutine [Humlicek,J Quant Spectrosc Radiat Transf 21 309 (1979)]
// for the calculation of the complex probability function
// ============================================================================================================================================================================================================================ //

long simulation::doHumlicek(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return doLorentz(area,y,ptrParam);

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if (!lHWHM)
		return doDoppler(area,y,ptrParam);

	std::vector<double> paramProfile(2,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;
	const double center = *(ptrParam+Wavelength) + delta;
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = fabs(lHWHM*invGWidth);

	const double amp = (!area) ? (amplitude/humlicekProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*humlicekProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//****************************************//
																						// VOIGT PROFILE (HUMLICEK APPROXIMATION) //
																						//****************************************//

// ============================================================================================================================================================================================================================ //

double simulation::humlicekProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	dcplx z(x,(*(ptrParamP+1)));
	dcplx cwerf = CPF(z) ;
	double intensity = (*ptrParamP)*cwerf.real();

    return intensity;
}

// ============================================================================================================================================================================================================================ //

																						//*********************************************//
																						// AREA VOIGT PROFILE (HUMLICEK APPROXIMATION) //
																						//*********************************************//
 
// ============================================================================================================================================================================================================================ //

double simulation::doAreaHumlicek(std::vector<double>::const_iterator ptrParam)
{
	return doAreaVoigt(ptrParam);
}

// ============================================================================================================================================================================================================================ //
//
// Hard collision model for velocity-changing (noted VC) collisions:
// In this case, the velocity memory is lost after each collision.
// In this so-called hard (labeled H) collision approximation,
//   H
// nu    (noted nuH_VC) is the hard collision rate, which is speed-independent.
//   VC
//  It is equal to the kinetic rate nuKin since each VC collision is thus assumed
//  to randomize the radiator velocity. Let us mention that the hard collision
//  approximation does not imply m p /m >> 1, as the simple model of rigid spheres,
//  and the usual picture of light radiators and heavy perturbers is somewhat excessive
//  for hard collisions.
//
//  In dimensionless variable: eta = nuH_VC/DWd = nuKin/DWd
//
// (*) In "true" gas media, the optically active molecule (the "radiator") cannot be con-
// sidered as alone, so that knowledge of the previously introduced intrinsic spectroscopic
// parameters is not sufficient for the modeling of the spectra or time-dependent signals.
// Except  for  very  low  pressures  where  considering  only  the  Doppler  effect  may  be
// sufficient, one must take into account the fact that the radiator is diluted in a "bath" of
// molecules or atoms (the "perturbers") with which it can interact
//
// ============================================================================================================================================================================================================================ //

																								//**********************************//
																								// NELKIN-GHATAK PROFILE SIMULATION //
																								//**********************************//

// ============================================================================================================================================================================================================================ //

long simulation::doNelkinGhatak(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));

	if (!lHWHM && !VCF)
		return doDoppler(area,y,ptrParam);
	if (!VCF)
		return doVoigt(area,y,ptrParam);

	std::vector<double> paramProfile(3,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;
	const double center = *(ptrParam+Wavelength) + delta;
	const double eta = fabs(VCF*invGWidth);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = lHWHM*invGWidth + eta;
	paramProfile[2] = SQRTPI*eta;

	const double amp = (!area) ? (amplitude/nelkinGhatakProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*nelkinGhatakProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																						//***********************//
																						// NELKIN-GHATAK PROFILE //
																						//***********************//

// ============================================================================================================================================================================================================================ //

double simulation::nelkinGhatakProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	dcplx z(x,(*(ptrParamP+1)));
	dcplx cwerf = complexErrFunc(z); // cf. voigt.cpp pour la routine "complexErrFunc"
	dcplx ratio = (cwerf*fcnInverse((double)1.0-(*(ptrParamP+2))*cwerf));

	return (*(ptrParamP)*ratio.real());
}

// ============================================================================================================================================================================================================================ //

																						//****************************//
																						// AREA NELKIN-GHATAK PROFILE //
																						//****************************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaNelkinGhatak(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return dblERREUR;
	if(!amplitude)
		return dblOK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));

	if (!lHWHM && !VCF)
		return doAreaDoppler(ptrParam);
	if (!VCF)
		return doAreaVoigt(ptrParam);

	std::vector<double> paramProfile(3,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;
	const double center = *(ptrParam+Wavelength) + delta;
	const double eta = fabs(VCF*invGWidth);
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = lHWHM*invGWidth + eta;
	paramProfile[2] = SQRTPI*eta;

	double norm = fcnInverse(nelkinGhatakProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //

																						//***********************************************************//
																						// NELKIN-GHATAK PROFILE SIMULATION (HUMLICEK APPROXIMATION) //
																						//***********************************************************//

// ============================================================================================================================================================================================================================ //

long simulation::doNelkinGhatakHumlicek(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));

	if (!lHWHM && !VCF)
		return doDoppler(area,y,ptrParam);
	if (!VCF)
		return doVoigt(area,y,ptrParam);

	std::vector<double> paramProfile(3,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;
	const double center = *(ptrParam+Wavelength) + delta;
	const double eta = fabs(VCF*invGWidth);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = lHWHM*invGWidth + eta;
	paramProfile[2] = SQRTPI*eta;

	const double amp = (!area) ? (amplitude/nelkinGhatakHumlicekProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*nelkinGhatakHumlicekProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																						//************************************************//
																						// NELKIN-GHATAK PROFILE (HUMLICEK APPROXIMATION) //
																						//************************************************//

// ============================================================================================================================================================================================================================ //

double simulation::nelkinGhatakHumlicekProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	dcplx z(x,(*(ptrParamP+1)));
	dcplx cwerf = CPF(z); // cf. humlicek.cpp pour la routine "CPF"
	dcplx ratio = (cwerf*fcnInverse((double)1.0-(*(ptrParamP+2))*cwerf));

	return (*(ptrParamP)*ratio.real());
}

// ============================================================================================================================================================================================================================ //

																						//*****************************************************//
																						// AREA NELKIN-GHATAK PROFILE (HUMLICEK APPROXIMATION) //
																						//*****************************************************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaNelkinGhatakHumlicek(std::vector<double>::const_iterator ptrParam)
{
	return doAreaNelkinGhatak(ptrParam);
}

// ============================================================================================================================================================================================================================ //

																						//****************************//
																						// GALATRY PROFILE SIMULATION //
																						//****************************//

// ============================================================================================================================================================================================================================ //

long simulation::doGalatry(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double invGWidth = fcnInverse(gHWHM);
	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;

	if ((!lHWHM) && (!VCF) && (!delta))
		return doDoppler(area,y,ptrParam);
	if (!VCF)
		return doVoigt(area,y,ptrParam);

	std::vector<double> paramProfile(5,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double center = *(ptrParam+Wavelength) + delta;
	const double yg = lHWHM*invGWidth;
	const double z = VCF*invGWidth;
	const double invz = fcnInverse(z);
	const double zFactor_X0 = (double)0.5*invz + yg;
	const double zzHGC = (double)0.5*fcnCarre(invz);
	const double bzHGC_X0 = (double)1.+ zFactor_X0*invz;
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_PI*invGWidth;
	paramProfile[1] = bzHGC_X0;
	paramProfile[2] = zFactor_X0;
	paramProfile[3] = zzHGC;
	paramProfile[4] = invz;

	const double amp = (!area) ? (amplitude/galatryProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*galatryProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//*****************//
																						// GALATRY PROFILE //
																						//*****************//

// ============================================================================================================================================================================================================================ //

double simulation::galatryProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	static const int lnchf = 0; //Log-Ausgabe (1)
    static const int ip = 10;   //Anzahl der genauen Stellen
	const double zFactor = *(ptrParamP+2);
	dcplx galatry((double)0.,(double)0.);
	dcplx denominateurCplxFactor(zFactor,-x);
	dcplx factor = fcnInverse(denominateurCplxFactor);
	MathFonctions::HypergeometriqueConfluente::complex rzHGC = {(double)0.,(double)0.};
	MathFonctions::HypergeometriqueConfluente::complex azHGC = {(double)1.0,(double)0.};
    MathFonctions::HypergeometriqueConfluente::complex bzHGC = {(*(ptrParamP+1)),-x*(*(ptrParamP+4))};
    MathFonctions::HypergeometriqueConfluente::complex zzHGC = {(*(ptrParamP+3)),(double)0.};
	MathFonctions::HypergeometriqueConfluente::conhyp_(&rzHGC, &azHGC, &bzHGC, &zzHGC, &lnchf, &ip);
	galatry = factor*dcplx(rzHGC.r,rzHGC.i);

	return ((*(ptrParamP))*galatry.real());

}

// ============================================================================================================================================================================================================================ //

																						//**********************//
																						// AREA GALATRY PROFILE //
																						//**********************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaGalatry(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return dblERREUR;
	if(!amplitude)
		return dblOK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));

	if ((!lHWHM) && (!VCF))
		return doAreaDoppler(ptrParam);
	if (!VCF)
		return doAreaVoigt(ptrParam);

	std::vector<double> paramProfile(5,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double yg = lHWHM*invGWidth;
	const double z = VCF*invGWidth;
	const double invz = fcnInverse(z);
	const double zFactor_X0 = (double)0.5*invz + yg;
	const double zzHGC = (double)0.5*fcnCarre(invz);
	const double bzHGC_X0 = (double)1.+ zFactor_X0*invz;
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_PI*invGWidth;
	paramProfile[1] = bzHGC_X0;
	paramProfile[2] = zFactor_X0;
	paramProfile[3] = zzHGC;
	paramProfile[4] = invz;

	double norm = fcnInverse(galatryProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //

																						//*************************************//
																						// INCOMPLETE GAMMA PROFILE SIMULATION //
																						//*************************************//

// ============================================================================================================================================================================================================================ //

long simulation::doIncompleteGamma(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double invGWidth = fcnInverse(gHWHM);
	const double VCF = fabs(*(ptrParam+velocityChanging));
	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	const double deltaP = (*(ptrParam+parameter_0));
	const double delta = (deltaP==deltaP) ? (deltaP*invGWidth):(double)0.0;

	if ((!lHWHM) && (!VCF) && (!delta))
		return doDoppler(area,y,ptrParam);
	if (!VCF)
		return doVoigt(area,y,ptrParam);

	std::vector<double>::const_iterator ptrFxRegion(fX.begin());
	std::vector<double>::const_iterator ptrFxEnd(fX.end());
	const double center = *(ptrParam+Wavelength) + delta;
	const double yg = lHWHM*invGWidth;
	const double z = VCF*invGWidth;
	double lineCenter((double)0.0);
	double *ptrY = &y[0];
	long region(0);

	// Region 1
	if ((z<(double)0.04) && (yg<(double)0.5))
    {
		//if all points are in region 1 set region = 1 else region = 0
		double br = (double)2.0*fabs(gHWHM);
		for(;ptrFxRegion!=ptrFxEnd;++ptrFxRegion)
			if (fabs((*ptrFxRegion)- center)>br)
					break;
		region = (ptrFxRegion==ptrFxEnd) ? 1 : 0;
	}
	//Region 2
    else if(((z>(double)0.1) && (yg<(double)4.0*pow(z,0.868))) || z>(double)5.)
		region = 2;
	//Region 3
    else
		region = 3;
	switch(region)
    {
		case 0:
				{
					std::vector<double>::const_iterator ptrFx(fX.begin());
					std::vector<double> paramProfile(15,(double)0.0);
					std::vector<double>::iterator ptrParamProfile = paramProfile.begin();
					long n = 2 + long((double)37.0*exp(-(double)0.6*yg));
					double amplitudeX(0.0);
					paramProfile[0] = INVERSE_PI*invGWidth;
					paramProfile[1] = yg;
					paramProfile[2] = z;
					paramProfile[3] = (double)n;
					paramProfile[4] = yg + paramProfile[3]*z;
					paramProfile[5] = (double)0.5*paramProfile[3];
					coeffAsymtoticExpansion(ptrParamProfile);
					const double amp = (!area) ? (amplitude/incompleteGammaRegionI(ptrParamProfile,lineCenter)):amplitude;
					for(;ptrFx!=ptrFxEnd;++ptrFx)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;
						amplitudeX = (lineCenter>(double)2.) ? incompleteGammaRegionIII(ptrParamProfile,lineCenter):incompleteGammaRegionI(ptrParamProfile,lineCenter);
						(*ptrY++)+=amp*amplitudeX;
					}
				}
				break;
		case 1:
				{
					std::vector<double>::const_iterator ptrFx(fX.begin());
					std::vector<double> paramProfile(15,(double)0.0);
					std::vector<double>::iterator ptrParamProfile = paramProfile.begin();
					paramProfile[0] = INVERSE_PI*invGWidth;
					paramProfile[1] = yg;
					paramProfile[2] = z;
					coeffAsymtoticExpansion(ptrParamProfile);
					const double amp = (!area) ? (amplitude/incompleteGammaRegionI(ptrParamProfile,lineCenter)):amplitude;
					for(;ptrFx!=ptrFxEnd;++ptrFx)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;						
						(*ptrY++)+=amp*incompleteGammaRegionI(ptrParamProfile,lineCenter);
					}
				}
				break;
		case 2:
				{
					std::vector<double>::const_iterator ptrFx(fX.begin());
					std::vector<double> paramProfile(6,(double)0.0);
					std::vector<double>::iterator ptrParamProfile = paramProfile.begin();
					long n = 4 + long(pow(z,-1.05)*((double)1.0+(double)3.0*exp(-(double)1.1*yg)));;
					paramProfile[0] = INVERSE_PI*invGWidth;
					paramProfile[1] = yg;
					paramProfile[2] = z;
					paramProfile[3] = (double)n;
					paramProfile[4] = (double)0.5/fcnCarre(z);
					paramProfile[5] = paramProfile[4] + yg/z;
					const double amp = (!area) ? (amplitude/incompleteGammaRegionII(ptrParamProfile,lineCenter)):amplitude;
					for(;ptrFx!=ptrFxEnd;++ptrFx)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;
						(*ptrY++)+=amp*incompleteGammaRegionII(ptrParamProfile,lineCenter);
					}
				}
				break;
		case 3:
				{
					std::vector<double>::const_iterator ptrFx(fX.begin());
					std::vector<double> paramProfile(6,(double)0.0);
					std::vector<double>::iterator ptrParamProfile = paramProfile.begin();
					long n = 2 + long((double)37.0*exp(-((double)0.6*yg)));
					paramProfile[0] = INVERSE_PI*invGWidth;
					paramProfile[1] = yg;
					paramProfile[2] = z;
					paramProfile[3] = (double)n;
					paramProfile[4] = yg + paramProfile[3]*z;
					paramProfile[5] = (double)0.5*paramProfile[3];
					const double amp = (!area) ? (amplitude/incompleteGammaRegionIII(ptrParamProfile,lineCenter)):amplitude;
					for(;ptrFx!=ptrFxEnd;++ptrFx)
					{
						lineCenter = ((*ptrFx)-center)*invGWidth;
						(*ptrY++)+=amp*incompleteGammaRegionIII(ptrParamProfile,lineCenter);
					}
				}
				break;
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																								//**************************//
																								// INCOMPLETE GAMMA PROFILE //
																								//**************************//
// [P.L. Varghese, R.K. Hanson, Collisional narrowing effects on spectral line shapes measured at high resolution.  Appl Opt 23, 2376 (1984)]
// ============================================================================================================================================================================================================================ //

// ============================================================================================================================================================== //

			// ************************************************************************************************************************************* //
			/*                         Calcul de l'amplitude de la fonction GALATRY dans la région I                                                 */
			/*                                               z<0.04 & yg<0.5                                                                         */
			// ************************************************************************************************************************************* //

// ============================================================================================================================================================== //

double simulation::incompleteGammaRegionI(std::vector<double>::iterator ptrCoefficientGalatry,const double &x)
{
	std::vector<double> wr(9,(double)0.0);
	std::vector<double> wi(9,(double)0.0);
	double yg(*(ptrCoefficientGalatry+1));
	double br = (double)2.0*x;
	double bi = (double)2.0*yg;
	dcplx z(x,yg);
	dcplx w = complexErrFunc(z);
	double valueW(0.0);

	wr[0] = w.real();
	wi[0] = w.imag();
	wr[1] = -br*wr[0] + bi*wi[0];
	wi[1] = TWObySQRTPI - br*wi[0]-bi*wr[0];
	for(long j=2,nd=2;j<9;++j,nd+=2)
	{
		wr[j] = -br*wr[j-1] + bi*wi[j-1] - (double)nd*wr[j-2];
		wi[j] = -br*wi[j-1] - bi*wr[j-1] - (double)nd*wi[j-2];
	}
	valueW = wr[0] - (*(ptrCoefficientGalatry+9))*wi[3] + (*(ptrCoefficientGalatry+10))*wr[4] + (*(ptrCoefficientGalatry+11))*wi[5]
	               - (*(ptrCoefficientGalatry+12))*wr[6] - (*(ptrCoefficientGalatry+13))*wi[7] + (*(ptrCoefficientGalatry+14))*wr[8];

	return (*(ptrCoefficientGalatry)*SQRTPI*valueW);
}

// ============================================================================================================================================================== //

			// ************************************************************************************************************************************* //
			/*                         Calcul de l'amplitude de la fonction GALATRY dans la région II                                                */
			/*                                    [z>0.1 & yg<4.0*pow(z,0.868)] ou z>5                                                               */                              
			// ************************************************************************************************************************************* //

// ============================================================================================================================================================== //

double simulation::incompleteGammaRegionII(std::vector<double>::iterator ptrCoefficientGalatry,const double &x)
{
	long n = (long)(*(ptrCoefficientGalatry+3));
	double invZ = fcnInverse(*(ptrCoefficientGalatry+2));
	double delta  = *(ptrCoefficientGalatry+4);
	double thetar = *(ptrCoefficientGalatry+5);
	double thetai = -x*invZ;
	double cr =  thetar/(fcnPythagore(thetar,thetai));
	double ci = -thetai/(fcnPythagore(thetar,thetai));
	double ct(0.0);
	double ar(cr), ay(0.0);
	double br(0.0),bi(0.0);

	for(long j=1;j<=n;++j)
	{
		ay = delta/(fcnPythagore(thetar+(double)j,thetai));
		br =  ay*(thetar+(double)j);
		bi = -ay*thetai;
		ct = br*cr-bi*ci;
		ci = bi*cr+br*ci;
		cr = ct;
		ar+=cr;
	}
	return ((*(ptrCoefficientGalatry)*ar)*invZ);
}

// ============================================================================================================================================================== //

			// ************************************************************************************************************************************* //
			/*                         Calcul de l'amplitude de la fonction GALATRY dans la région III                                               */
			/*                                          les cas restants...																			 */
			// ************************************************************************************************************************************* //

// ============================================================================================================================================================== //

double simulation::incompleteGammaRegionIII(std::vector<double>::iterator ptrCoefficientGalatry,const double &x)
{

	long n    = (long)(*(ptrCoefficientGalatry+3));
	double yg = *(ptrCoefficientGalatry+1);
	double z  = *(ptrCoefficientGalatry+2);
	double br = *(ptrCoefficientGalatry+4);
	double bi = -x;
	double nx = *(ptrCoefficientGalatry+5);
	double ax = nx/(fcnPythagore(br,bi));
	double ar =  br*ax;
	double ai = -bi*ax;
	for(long j=n-1;j>0;--j)
	{
		br-=z;
		nx-=(double)0.5;
		ax = nx/(fcnPythagore(br+ar,bi+ai));
		ar =  (ar+br)*ax;
		ai = -(ai+bi)*ax;
	}
	return  ((*(ptrCoefficientGalatry))*(yg+ar)*fcnInverse(fcnPythagore(yg+ar,bi+ai)));
}

// ============================================================================================================================================================================================================================ //

																						//*******************************//
																						// AREA INCOMPLETE GAMMA PROFILE //
																						//*******************************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaIncompleteGamma(std::vector<double>::const_iterator ptrParam)
{
	return doAreaGalatry(ptrParam);
}

// ============================================================================================================================================================================================================================ //


// ============================================================================================================================================================================================================================ //

																						//*************************//
																						// FANO PROFILE SIMULATION //
																						//*************************//

// ============================================================================================================================================================================================================================ //
																				/* According to U. Fano */
// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | Fano                                        				     	                |         asymmetry Factor q      |	                           |	                        |	                         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

long simulation::doFano(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(amplitude<(double)0.0)
		return ERREUR;
	if(!amplitude)
		return OK;

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if(!lHWHM)
		return OK;

	const double q = *(ptrParam+parameter_0);
	if (!q)
		return doLorentz(area,y,ptrParam);
	if (fabs(q)>=lHWHM)
		return ERREUR;

	std::vector<double> paramProfile(3,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double center = *(ptrParam+Wavelength);
	const double invlHWHM = fcnInverse(lHWHM);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_PI*invlHWHM*fcnInverse(1-fcnCarre(q*invlHWHM));
	paramProfile[1] = lHWHM;
	paramProfile[2] = q;

	const double amp = (!area) ? (amplitude/fanoProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
			lineCenter = ((*ptrFx)-center);
			(*ptrY++)+=amp*fanoProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																						//**************//
																						// FANO PROFILE //
																						//**************//

// ============================================================================================================================================================================================================================ //

double simulation::fanoProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	double yg = (*(ptrParamP+1));
	double lineCenterQ2 = fcnCarre(x+(*(ptrParamP+2)));
	double intensity = *(ptrParamP)*((double)1.0-lineCenterQ2*fcnInverse(fcnPythagore(x,yg)));

	return intensity;
}

// ============================================================================================================================================================================================================================ //

																						//*******************//
																						// AREA FANO PROFILE //
																						//*******************//

// ============================================================================================================================================================================================================================ //

double simulation::doAreaFano(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double lHWHM = fabs(*(ptrParam+gammaHWHM));
	if(!lHWHM)
		return OK;

	const double q = *(ptrParam+parameter_0);
	if (!q)
		return doAreaLorentz(ptrParam);

	std::vector<double> paramProfile(3,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double center = *(ptrParam+Wavelength);
	const double invlHWHM = fcnInverse(lHWHM);
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_PI*invlHWHM*fcnInverse(1-fcnCarre(q*invlHWHM));
	paramProfile[1] = lHWHM;
	paramProfile[2] = q;

	double norm = fanoProfile(ptrParamProfile,lineCenter);

	return (amplitude*norm);
}

// ============================================================================================================================================================================================================================ //

																						/*************************************/
																						/* SPEED-DEPENDENT LINE-SHAPE MODELS */
																						/*************************************/

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																						//******************************************************//
																						// a (useful) quadratic model for Gamma(v) and Delta(v) // 
																						//******************************************************//

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //
// In terms of dimensionless variables, the complex line shape function for
// the the speed-dependent Voigt profile (SDVP) is:
// Isdvp(u) = (2/PI^3/2)*$$\int_(-INF^+INF) exp(-x*x)*x*{arctg[f(u,x)])
//                                         +i/2*ln[1+f(u,x)*f(u,x)]}dx$$ (1)
// where u = (nu - nu0)/nuD is the reduced spectral detuning. The
// variable of integration, x = v/vp is the reduced absorber speed,
// vp = sqrt(2*kB*T/mA) , mA is the absorber mass, T is the gas  temperature,
// and  kB is  the  Boltzmann  constant.
// The Doppler width [full width at half maximum (FWHM)] is
// GammaD =2*sqrt(ln(2))nuD , and the Lorentzian (collisional) width [half
// width at half maximum (HWHM)] and shift of line center are
// Gamma and Delta, respectively. The function in the integrand of Eq. (1)
// is f(u,x) = [u - d*Bs (x) + x)]/[g*Bw(x)], where g = Gamma/nuD ,
// d = Delta/nuD , and Bw(x) and Bs(x)  are  the  reduced  speed-dependent collisional
// width and shift functions, respectively.
// Assuming  a  quadratic  dependence  for  each  of  the  speed-dependent
// collisional functions then
//                Bw(x) = 1 + aw*(x^2 - 3/2),
//				  Bs(x) = 1 + as*(x^2 - 3/2).
// define the broadening speed-dependence parameter aw, and
// the shifting speed-dependence parameter as. These two
// terms are treated as fitted parameters in the data analysis. At
// fixed temperature, aw and as are expected to be independent
// of gas pressure.
//
// ============================================================================================================================================================================================================================ //
//
// Reference:
// J Chem Phys. 2011 Aug 14;135(6):064308. doi: 10.1063/1.3624527.
// The air-broadened, near-infrared CO2 line shape in the spectrally isolated regime: evidence of simultaneous Dicke narrowing and speed dependence.
// Long DA, Bielska K, Lisak D, Havey DK, Okumura M, Miller CE, Hodges JT.
//
// ============================================================================================================================================================================================================================ //

																						//****************************************************//
																						// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE SIMULATION //
																						//****************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3				 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed Dependent Voigt  	   					                            |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentVoigt(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double aw = fabs(*(ptrParam+parameter_1));
	const double as = *(ptrParam+parameter_2);

	if((!aw) && (!as))
		return doVoigt(area,y,ptrParam);

	const double delta = *(ptrParam+parameter_0);
	const double gamma = fabs(*(ptrParam+gammaHWHM));
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const long SDVProfile(5); 
	const long nbParam(SDVProfile+3*NB_PTS_GAUSS_LEGENDRE);
	double *ptrY = &y[0];
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = TWObyCUBERTPI*invGWidth;
	paramProfile[1] = gamma*invGWidth;
	paramProfile[2] = delta*invGWidth;
	paramProfile[3] = aw;
	paramProfile[4] = as;
	initQSDVP(SDVProfile,ptrParamProfile); 
	
	const double amp = (!area) ? (amplitude/quadraticSpeedDependentVoigtProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*quadraticSpeedDependentVoigtProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																						//*****************************************//
																						// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE //
																						//*****************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3				 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed Dependent Voigt  	   					                            |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentVoigtProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	const long shiftSDV = 5;
	std::vector<double>::const_iterator ptrBs(ptrParamP+shiftSDV);
	std::vector<double>::const_iterator ptrBw(ptrBs+NB_PTS_GAUSS_LEGENDRE);
	std::vector<double>::const_iterator ptrIntFcn(ptrBw+NB_PTS_GAUSS_LEGENDRE);
	double integraleSDV_Pos((double)0.0);
	double integraleSDV_Neg((double)0.0);
	double starPtrBs(0.0);
	double starPtrBw(0.0);
	double starPtrIntFcn(0.0);
	double xGL(0.0);	

	for(long index=0;index<NB_PTS_GAUSS_LEGENDRE;++index)
	{
		xGL = PIby2*x1024[index];
		starPtrBs = (*ptrBs++);
		starPtrBw = (*ptrBw++);
		starPtrIntFcn = (*ptrIntFcn++);
		integraleSDV_Pos +=  starPtrIntFcn*atan((tan(xGL)-starPtrBs+x)/starPtrBw);
		integraleSDV_Neg += -starPtrIntFcn*atan((tan(-xGL)-starPtrBs+x)/starPtrBw);		
	}

	return (PIby2*(*(ptrParamP))*(integraleSDV_Pos+integraleSDV_Neg));
}

// ============================================================================================================================================================================================================================ //

																						//**********************************************//
																						// AREA QUADRATIC SPEED-DEPENDENT VOIGT PROFILE //
																						//**********************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3				 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed Dependent Voigt  	   					                            |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentVoigt(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double aw = fabs(*(ptrParam+parameter_1));
	const double as = *(ptrParam+parameter_2);

	if((!aw) && (!as))
		return doAreaVoigt(ptrParam);

	const double delta = *(ptrParam+parameter_0);
	const double gamma = fabs(*(ptrParam+gammaHWHM));
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const long SDVProfile(5); 
	const long nbParam(SDVProfile+3*NB_PTS_GAUSS_LEGENDRE);
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = TWObyCUBERTPI*invGWidth;
	paramProfile[1] = gamma*invGWidth;
	paramProfile[2] = delta*invGWidth;
	paramProfile[3] = aw;
	paramProfile[4] = as;
	initQSDVP(SDVProfile,ptrParamProfile);

	const double norm = fcnInverse(quadraticSpeedDependentVoigtProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// The complex speed-dependent Nelkin-Ghatak profile (SDNGP) can be expressed in
// terms of the SDVP by
// Isdngp(u) = Isdvp*(u)/(1-PI*z*Isdvp*(u))
// where  z = nuNAR /nuD and  Isdvp*(u)
// is  the  function  Isdvp(u) with gBw(x) + z
// substituted for gBw(x). Here, the narrowing
// frequency, nuNAR, accounts for the collisional narrowing effect,
// and at fixed temperature it is proportional to pressure.
//
// ============================================================================================================================================================================================================================ //
//
// Reference:
// J Chem Phys. 2011 Aug 14;135(6):064308. doi: 10.1063/1.3624527.
// The air-broadened, near-infrared CO2 line shape in the spectrally isolated regime: evidence of simultaneous Dicke narrowing and speed dependence.
// Long DA, Bielska K, Lisak D, Havey DK, Okumura M, Miller CE, Hodges JT.
//
// ============================================================================================================================================================================================================================ //

																				//************************************************************//
																				// QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE SIMULATION //
																				//************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Nelkin-Ghatak   					                        |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentNelkinGhatak(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double dzeta = fabs(*(ptrParam+velocityChanging));
	if(!dzeta)
		return doQuadraticSpeedDependentVoigt(area,y,ptrParam);

	const double aw = fabs(*(ptrParam+parameter_1));
	const double as = *(ptrParam+parameter_2);

	if((!aw) && (!as))
		return doNelkinGhatak(area,y,ptrParam);

	const double delta = *(ptrParam+parameter_0);
	const double gamma = fabs(*(ptrParam+gammaHWHM));
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const long SDNGProfile(6); 
	const long nbParam(SDNGProfile+3*NB_PTS_GAUSS_LEGENDRE);
	double *ptrY = &y[0];
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = TWObyCUBERTPI*invGWidth;
	paramProfile[1] = gamma*invGWidth;
	paramProfile[2] = delta*invGWidth;
	paramProfile[3] = aw;
	paramProfile[4] = as;
	paramProfile[5] = dzeta*invGWidth;
	initQSDNGP(SDNGProfile,ptrParamProfile);

	const double amp = (!area) ? (amplitude/quadraticSpeedDependentNelkinGhatakProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*quadraticSpeedDependentNelkinGhatakProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																				//*************************************************//
																				// QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE //
																				//*************************************************//

// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentNelkinGhatakProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	const long shiftSDNG = 6;
	const double zFactor = TWObySQRTPI*(*(ptrParamP+5)); // PI*z*Isdvp*(u) avec Isdvp*(u)= (2/PI^3/2)*$$\int_(-INF^+INF) exp(-x*x)*x*{arctg[f(u,x)])+i/2*ln[1+f(u,x)*f(u,x)]}dx$$
	double integraleSDV_Real((double)0.0);
	double integraleSDV_Imag((double)0.0);
	double argPos((double)0.0);
	double argNeg((double)0.0);
	double xGL(0.0);
	double starPtrBs(0.0);
	double starPtrBw(0.0);
	double starPtrIntFcn(0.0);
	std::vector<double>::const_iterator ptrBs(ptrParamP+shiftSDNG);
	std::vector<double>::const_iterator ptrBw(ptrBs+NB_PTS_GAUSS_LEGENDRE);
	std::vector<double>::const_iterator ptrIntFcn(ptrBw+NB_PTS_GAUSS_LEGENDRE);
	dcplx ISDVStar((double)0.0,(double)0.0);
	dcplx ISDNG((double)0.0,(double)0.0);

	for(long index=0;index<NB_PTS_GAUSS_LEGENDRE;++index)
	{
		xGL = PIby2*x1024[index];
		starPtrBs = (*ptrBs++);
		starPtrBw = (*ptrBw++);
		starPtrIntFcn = (*ptrIntFcn++);
		argPos = (x-starPtrBs+tan(xGL))/starPtrBw;
		argNeg = (x-starPtrBs-tan(xGL))/starPtrBw;
		integraleSDV_Real += starPtrIntFcn*(atan(argPos)-atan(argNeg));
		integraleSDV_Imag += starPtrIntFcn*(double)0.5*(log((double)1.+fcnCarre(argPos))-log((double)1.+fcnCarre(argNeg)));
	}

	ISDVStar = PIby2*dcplx(integraleSDV_Real,integraleSDV_Imag);
	ISDNG = (*(ptrParamP)*ISDVStar)/((double)1.0-zFactor*ISDVStar);

	return (ISDNG.real());

}

// ============================================================================================================================================================================================================================ //

																				//******************************************************//
																				// AREA QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE //
																				//******************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Nelkin-Ghatak   					                        |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentNelkinGhatak(std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double dzeta = *(ptrParam+velocityChanging);
	if(!dzeta)
		return doAreaQuadraticSpeedDependentVoigt(ptrParam);

	const double aw = fabs(*(ptrParam+parameter_1));
	const double as = *(ptrParam+parameter_2);

	if((!aw) && (!as))
		return doAreaNelkinGhatak(ptrParam);

	const double delta = *(ptrParam+parameter_0);
	const double gamma = fabs(*(ptrParam+gammaHWHM));
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const long SDNGProfile(6); 
	const long nbParam(SDNGProfile+3*NB_PTS_GAUSS_LEGENDRE);
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = TWObyCUBERTPI*invGWidth;
	paramProfile[1] = gamma*invGWidth;
	paramProfile[2] = delta*invGWidth;
	paramProfile[3] = aw;
	paramProfile[4] = as;
	paramProfile[5] = fabs(dzeta)*invGWidth;
	initQSDNGP(SDNGProfile,ptrParamProfile);

	const double norm = fcnInverse(quadraticSpeedDependentNelkinGhatakProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// In  Ref[2],  it  was  shown  that  the  quadratic-speed-dependent  Voigt  profile  can  be
// expressed as a combination of two "usual" Voigt functions.
// The quadratic-speed-dependent (complex) Voigt line profile becomes:
// Isdvp(u) = c/(sqrt(PI)*u0*vp)*[w(iZ1)-w(iZ2)]
// with
// Z1 = sqrt(X+Y) - sqrt(Y)
// Z2 = sqrt(X+Y) + sqrt(Y)
// vp = sqrt(2*kB*T/m), the most probable speed
// X = [i(u-u0)+C0]/C2
// Y = (u0*vp/(2*c*C2))^2 = (gammaD/(2*C2*sqrt(ln2)))^2 (gammaD is the Doppler width)
// gamma(v)-i*delta(v) = C0 + C2*{(v/vp)^2-1.5}
// C0 = gamma0 - i*delta0
// C2 = (gamma2 - i*delta2)
//
// ============================================================================================================================================================================================================================ //
// References:
// [1] NH Ngo, D Lisak, H Tran, JM Hartmann. An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes
//     J Quant Spectrosc Radiat Transf 2013
// [2] Boone CD, Walker KA, Bernath PF. Speed-dependent Voigt profile for water vapor in infrared remote sensing applications.
//     J Quant Spectrosc Radiat Transf 2007;105:525-532.
//
// ============================================================================================================================================================================================================================ //

																		//******************************************************************//
																		// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE SIMULATION (BOONE MODEL) //
																		//******************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Voigt (Boone model)		 		                        |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//  In this case, the  complex  probability function w(z) is computed thanks to the CERNlib routine
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentVoigtBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doVoigt(area,y,ptrParam);

	std::vector<double> paramProfile(9,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)*gHWHM;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();

	const double amp = (!area) ? (amplitude/quadraticSpeedDependentVoigtBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*quadraticSpeedDependentVoigtBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																				//*******************************************************//
																				// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE (BOONE MODEL) //
																				//*******************************************************//

// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentVoigtBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	static const dcplx iCplx(0.0,1.0);
	dcplx xCplx = (iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*dcplx(*(ptrParamP+7),*(ptrParamP+8));
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx iz1CPlx = iCplx*(sqrt(xCplx+dcplx(*(ptrParamP+5),*(ptrParamP+6)))-yCplx);
	dcplx iz2CPlx = iz1CPlx + iCplx*(double)2.0*yCplx;
	dcplx AtermByPI = (*ptrParamP)*(complexErrFunc(iz1CPlx)-complexErrFunc(iz2CPlx));

	return (AtermByPI.real());
}

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3				 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed Dependent Voigt (Boone model)			                            |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentVoigtBoone(std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doAreaVoigt(ptrParam);

	std::vector<double> paramProfile(9,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)*gHWHM;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();

	const double norm = fcnInverse(quadraticSpeedDependentVoigtBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// In  Ref[2],  it  was  shown  that  the  quadratic-speed-dependent  Voigt  profile  can  be
// expressed as a combination of two "usual" Voigt functions.
// The quadratic-speed-dependent (complex) Voigt line profile becomes:
// Isdvp(u) = c/(sqrt(PI)*u0*vp)*[w(iZ1)-w(iZ2)]
// with
// Z1 = sqrt(X+Y) - sqrt(Y)
// Z2 = sqrt(X+Y) + sqrt(Y)
// vp = sqrt(2*kB*T/m), the most probable speed
// X = [i(u-u0)+C0]/C2
// Y = (u0*vp/(2*c*C2))^2 = (gammaD/(2*C2*sqrt(ln2)))^2 (gammaD is the Doppler width)
// gamma(v)-i*delta(v) = C0 + C2*{(v/vp)^2-1.5}
// C0 = gamma0 - i*delta0
// C2 = (gamma2 - i*delta2)
//
// Within the Speed-Dependent Hard-Collision  model, the  spectral  profile  of an isolated
// line is given by can be written as:
//      Isdngp(u) = Isdvp(u)/(1-nuVC*PI*Isdvp(u))
// with nuVC the frequency of velocity-changing collisions when assuming no correlation
// between   velocity   and   rotational-state   changes
//
// ============================================================================================================================================================================================================================ //
// References:
// [1] NH Ngo, D Lisak, H Tran, JM Hartmann. An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes
//     J Quant Spectrosc Radiat Transf 2013
// [2] Boone CD, Walker KA, Bernath PF. Speed-dependent Voigt profile for water vapor in infrared remote sensing applications.
//     J Quant Spectrosc Radiat Transf 2007;105:525-532.
//
// ============================================================================================================================================================================================================================ //

																				//****************************************************************//
																				// QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (BOONE MODEL) //
																				//****************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Hard Collision (Boone model)				                |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//  In this case, the  complex  probability function w(z) is computed thanks to the CERNlib routine
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentNelkinGhatakBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double dzeta = *(ptrParam+velocityChanging);
	if(!dzeta)
		return doQuadraticSpeedDependentVoigtBoone(area,y,ptrParam);

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doNelkinGhatak(area,y,ptrParam);


	std::vector<double> paramProfile(10,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)/invGWidth;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = c0Tilde.real()+dzeta;
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();
	paramProfile[9] = dzeta;

	const double amp = (!area) ? (amplitude/quadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*quadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																				//****************************************************************//
																				// QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (BOONE MODEL) //
																				//****************************************************************//

// ============================================================================================================================================================================================================================ //
//	"qSDNGP_Humlicek": quadratic-Speed-Dependent Hard-Collision
//	Subroutine to Compute the complex normalized spectral shape of an
//	isolated line by the qSDHC mode
//
//	Called Routines: 'complexErrFunc'	(complex error function~Complex Probability Function)
// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentNelkinGhatakBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	static const dcplx iCplx(0.0,1.0);
	dcplx xCplx = (iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*dcplx(*(ptrParamP+7),*(ptrParamP+8));
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx iz1CPlx = iCplx*(sqrt(xCplx+dcplx(*(ptrParamP+5),*(ptrParamP+6)))-yCplx);
	dcplx iz2CPlx = iz1CPlx + iCplx*(double)2.0*yCplx;
	dcplx AtermByPI = (*ptrParamP)*(complexErrFunc(iz1CPlx)-complexErrFunc(iz2CPlx));
	dcplx cplxNG = AtermByPI/((double)1.0-(*(ptrParamP+9))*PI*AtermByPI);

	return (cplxNG.real());
}

// ============================================================================================================================================================================================================================ //

																				//**************************************************************************//
																				// AREA QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (HUMLICEK ROUTINE) //
																				//**************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Hard Collision (Boone model) 			                |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentNelkinGhatakBoone(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double dzeta = *(ptrParam+velocityChanging);
	if(!dzeta)
		return doAreaQuadraticSpeedDependentVoigtBoone(ptrParam);

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doAreaNelkinGhatak(ptrParam);

	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	std::vector<double> paramProfile(10,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)/invGWidth;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_SQRTPI*invGWidth;
	paramProfile[1] = c0Tilde.real()+dzeta;
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();
	paramProfile[9] = dzeta;

	const double norm = fcnInverse(quadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //

																				//*****************************************************************************//
																				// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//*****************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Voigt (Boone model with Humlicek routine)                |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //
// Based on a slightly improved version of the CPF subroutine [Humlicek,J Quant Spectrosc Radiat Transf 21 309 (1979)]
// for the calculation of the complex probability function
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentHumlicekBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doHumlicek(area,y,ptrParam);

	std::vector<double> paramProfile(9,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)*gHWHM;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = INVERSE_PI*invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();

	const double amp = (!area) ? (amplitude/quadraticSpeedDependentHumlicekBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*quadraticSpeedDependentHumlicekBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

																				//*****************************************************************************//
																				// QUADRATIC SPEED-DEPENDENT VOIGT PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//*****************************************************************************//

// ============================================================================================================================================================================================================================ //
//	"qSDVP_Humlicek": quadratic-Speed-Dependent Voigt
//	Subroutine to Compute the complex normalized spectral shape of an
//	isolated line by the qSDV model
//
//	Called Routines: 'CPF'	(Complex Probability Function)
//	---------------  'CPF3'	(Complex Probability Function for the region 3)
// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentHumlicekBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	static const dcplx iCplx((double)0.0,(double)1.0);
	dcplx invC2(*(ptrParamP+7),*(ptrParamP+8));
	dcplx xCplx = (iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*invC2;
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx y2Cplx(*(ptrParamP+5),*(ptrParamP+6));
	dcplx Aterm(1.0,1.0);
	dcplx z1CPlx(0.0,0.0);
	dcplx z2CPlx(0.0,0.0);
	dcplx w1CPlx(0.0,0.0);
	dcplx w2CPlx(0.0,0.0);
	double borneSup((double)3.0e-8*abs(y2Cplx));
	double center(0.0);
	double xz1(0.0);
	double xz2(0.0);
	double yz1(0.0);
	double yz2(0.0);
	double SZ1(0.0);
	double SZ2(0.0);
	double DSZ(0.0);
	double SZmax(0.0);
	double SZmin(0.0);
	long region(0);

	// when abs(X) is much larger than abs(Y)
	if (abs(y2Cplx)<(double)1.0e-15*abs(xCplx))
		region = 1;
	// when abs(Y) is much larger than abs(X)
	else if (abs(xCplx)<borneSup)
		region = 2;
	switch(region)
	{
		case 0:{
					// calculating Z1 and Z2
					z1CPlx = sqrt(xCplx+y2Cplx)-yCplx;
					z2CPlx = z1CPlx + (double)2.0*yCplx;
					// calculating the real and imaginary parts of Z1 and Z2
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					// check if Z1 and Z2 are close to each other
					SZ1 = fcnNorme(xz1,yz1);
					SZ2 = fcnNorme(xz2,yz2);
					DSZ = fabs(SZ1-SZ2);
					SZmax = max(SZ1,SZ2);
					SZmin = min(SZ1,SZ2);
					// when Z1 and Z2 are close to each other, ensure that they are in
					// the same interval of CPF
					if (DSZ<=(double)1.0 && SZmax>(double)8.0 &&  SZmin<=(double)8.0)
					{
						w1CPlx = CPF3(dcplx(xz1,yz1));
						w2CPlx = CPF3(dcplx(xz2,yz2));
					}
					else
					{
						w1CPlx = CPF(dcplx(xz1,yz1));
						w2CPlx = CPF(dcplx(xz2,yz2));
					}
					// calculating the A term of the profile
					Aterm = SQRTPI*(w1CPlx-w2CPlx);
			   }
			   break;
		case 1:{
					// when abs(X) is much larger than abs(Y)
					if (abs(sqrt(xCplx))<(double)4000.0)
					{
						xz1 = -(sqrt(xCplx)).imag();
						yz1 = (sqrt(xCplx)).real();
						w1CPlx = CPF(dcplx(xz1,yz1));
						Aterm = ((double)2.0*SQRTPI*invC2)*(INVERSE_SQRTPI-sqrt(xCplx)*w1CPlx);
					}
					else
					// when abs(X) is much larger than 1
					{
						Aterm = invC2*((double)1.0/xCplx-(double)1.0/fcnCarre(xCplx));
					}
			   }
			   break;
		case 2:{
					// when abs(Y) is much larger than abs(X)
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = (iCplx*x+c0t)*(*(ptrParamP));
					z2CPlx = sqrt(xCplx+y2Cplx) + sqrt(y2Cplx);
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					w1CPlx = CPF(dcplx(xz1,yz1));
					w2CPlx = CPF(dcplx(xz2,yz2));
					Aterm = SQRTPI*(w1CPlx-w2CPlx);
			   }
			   break;
		default:break;
	}

	return (*(ptrParamP)*Aterm.real());
}

// ============================================================================================================================================================================================================================ //

																				//***********************************************************************************//
																				// AREA QUADRATIC SPEED-DEPENDENT VOIGT PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//**********************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// quadratic Speed-Dependent Voigt (Boone model with Humlicek routine) 	                |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //
// Based on a slightly improved version of the CPF subroutine [Humlicek,J Quant Spectrosc Radiat Transf 21 309 (1979)]
// for the calculation of the complex probability function
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentHumlicekBoone(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doAreaHumlicek(ptrParam);

	std::vector<double> paramProfile(9,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = (double)1.0/c2;
	const dcplx yCplx = ((double)0.5*invC2)/invGWidth;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);

	paramProfile[0] = INVERSE_PI*invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();

	const double norm = fcnInverse(quadraticSpeedDependentHumlicekBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //

																				//**************************************************************************************//
																				// QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//**************************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Hard Collision (Boone model with Humlicek routine)       |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //
// Based on a slightly improved version of the CPF subroutine [Humlicek,J Quant Spectrosc Radiat Transf 21 309 (1979)]
// for the calculation of the complex probability function
// ============================================================================================================================================================================================================================ //

long simulation::doQuadraticSpeedDependentNelkinGhatakHumlicekBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double dzeta = fabs(*(ptrParam+velocityChanging));
	if(!dzeta)
		return doQuadraticSpeedDependentHumlicekBoone(area,y,ptrParam);

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doNelkinGhatakHumlicek(area,y,ptrParam);

	std::vector<double> paramProfile(10,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)*gHWHM;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real()+dzeta;
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();
	paramProfile[9] = dzeta;

	const double amp = (!area) ? (amplitude/quadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*quadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;

}

// ============================================================================================================================================================================================================================ //

																				//**************************************************************************************//
																				// QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//**************************************************************************************//

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //
//	"qSDNGP_Humlicek": quadratic-Speed-Dependent Hard-Collision
//	Subroutine to Compute the complex normalized spectral shape of an
//	isolated line by the qSDHC mode
//
//	Called Routines: 'CPF'	(Complex Probability Function)
//	---------------  'CPF3'	(Complex Probability Function for the region 3)
// ============================================================================================================================================================================================================================ //

double simulation::quadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	static const dcplx iCplx((double)0.0,(double)1.0);
	const double cte(*(ptrParamP));
	dcplx invC2(*(ptrParamP+7),*(ptrParamP+8));
	dcplx xCplx = (iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*invC2;
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx y2Cplx(*(ptrParamP+5),*(ptrParamP+6));
	dcplx Aterm(0.0,0.0);
	dcplx qSDNGP(0.0,0.0);
	dcplx z1CPlx(0.0,0.0);
	dcplx z2CPlx(0.0,0.0);
	dcplx w1CPlx(0.0,0.0);
	dcplx w2CPlx(0.0,0.0);
	double borneSup((double)3.0e-8*abs(y2Cplx));
	double center(0.0);
	double xz1(0.0);
	double xz2(0.0);
	double yz1(0.0);
	double yz2(0.0);
	double SZ1(0.0);
	double SZ2(0.0);
	double DSZ(0.0);
	double SZmax(0.0);
	double SZmin(0.0);
	long region(0);

	// when abs(X) is much larger than abs(Y)
	if (abs(y2Cplx)<(double)1.0e-15*abs(xCplx))
		region = 1;
	// when abs(Y) is much larger than abs(X)
	else if (abs(xCplx)<borneSup)
		region = 2;
	switch(region)
	{
		case 0:{
					// calculating Z1 and Z2
					z1CPlx = sqrt(xCplx+y2Cplx)-yCplx;
					z2CPlx = z1CPlx + (double)2.0*yCplx;
					// calculating the real and imaginary parts of Z1 and Z2
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					// check if Z1 and Z2 are close to each other
					SZ1 = fcnNorme(xz1,yz1);
					SZ2 = fcnNorme(xz2,yz2);
					DSZ = fabs(SZ1-SZ2);
					SZmax = max(SZ1,SZ2);
					SZmin = min(SZ1,SZ2);
					// when Z1 and Z2 are close to each other, ensure that they are in
					// the same interval of CPF
					if (DSZ<=(double)1.0 && SZmax>(double)8.0 &&  SZmin<=(double)8.0)
					{
						w1CPlx = CPF3(dcplx(xz1,yz1));
						w2CPlx = CPF3(dcplx(xz2,yz2));
					}
					else
					{
						w1CPlx = CPF(dcplx(xz1,yz1));
						w2CPlx = CPF(dcplx(xz2,yz2));
					}
					// calculating the A term of the profile
					Aterm = SQRTPI*cte*(w1CPlx-w2CPlx);
			   }
			   break;
		case 1:{
					// when abs(X) is much larger than abs(Y)
					if (abs(sqrt(xCplx))<(double)4000.0)
					{
						xz1 = -(sqrt(xCplx)).imag();
						yz1 = (sqrt(xCplx)).real();
						w1CPlx = CPF(dcplx(xz1,yz1));
						Aterm = ((double)2.0*SQRTPI*invC2)*(INVERSE_SQRTPI-sqrt(xCplx)*w1CPlx);
					}
					else
					// when abs(X) is much larger than 1
					{
						Aterm = invC2*((double)1.0/xCplx-(double)1.0/fcnCarre(xCplx));
					}
			   }
			   break;
		case 2:{
					// when abs(Y) is much larger than abs(X)
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = (iCplx*x+c0t)*(*(ptrParamP));
					z2CPlx = sqrt(xCplx+y2Cplx) + sqrt(y2Cplx);
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					w1CPlx = CPF(dcplx(xz1,yz1));
					w2CPlx = CPF(dcplx(xz2,yz2));
					Aterm = SQRTPI*cte*(w1CPlx-w2CPlx);
			   }
			   break;
		default:break;
	}

	qSDNGP = Aterm/((double)1.0-(*(ptrParamP+9))*Aterm);

	return (INVERSE_PI*qSDNGP.real());
}

// ============================================================================================================================================================================================================================ //

																				//**************************************************************************//
																				// AREA QUADRATIC SPEED-DEPENDENT HARD-COLLISION PROFILE (HUMLICEK ROUTINE) //
																				//**************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | quadratic Speed-Dependent Hard Collision (Boone model) 			                |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |			                 |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //
// Boone model: the  quadratic-speed-dependent Nelkin-Ghatak profile can be expressed as a combination of quadratic Speed Dependent Voigt (Boone model)
// ============================================================================================================================================================================================================================ //

double simulation::doAreaQuadraticSpeedDependentNelkinGhatakHumlicekBoone(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double dzeta = *(ptrParam+velocityChanging);
	if(!dzeta)
		return doAreaQuadraticSpeedDependentHumlicekBoone(ptrParam);

	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);

	if((!gamma2) && (!delta2))
		return doAreaNelkinGhatakHumlicek(ptrParam);

	std::vector<double> paramProfile(10,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double gamma0 = fabs(*(ptrParam+gammaHWHM));
	const double delta0 = *(ptrParam+parameter_0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = (c0 - (double)1.5*c2);
	const dcplx invC2 = fcnInverse(c2);
	const dcplx yCplx = ((double)0.5*invC2)/invGWidth;
	const dcplx y2Cplx = fcnCarre(yCplx);
	double lineCenter((double)0.0);

	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real()+dzeta;
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = invC2.real();
	paramProfile[8] = invC2.imag();
	paramProfile[9] = dzeta;

	const double norm = fcnInverse(quadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// In  Ref[2],  it  was  shown  that  the  quadratic-speed-dependent  Voigt  profile  can  be
// expressed as a combination of two "usual" Voigt functions.
// Within the  partially-Correlated Speed-Dependent Hard-Collision model of an isolated
// line is given by can be written as:
//      IpCqSDHC(u) = 1/PI*{A(u)/([1-nuVC-eta(C0-1.5*C2)]*A(u) + [eta*C2/vp^2)*B(u)])}
// The A(u) and B(u) terms can be expressed as combination of the complex probability function w(z):
// 1/PI*A(u) = c/(sqrt(PI)*u0*vp)*[w(iZ1)-w(iZ2)]
// B(u) = (vp^2/C2Tilde)*[-1+sqrt(PI)/(2*sqrt(Y)*{(1-Z1^2)*w(iZ1)-(1-Z2^2)*w(iZ2)})]
// with
// X = [i(u-u0)+C0Tilde]/C2Tilde
// Y = (u0*vp/(2*c*C2Tilde))^2 = (gammaD/(2*C2*sqrt(ln2)))^2 (gammaD is the Doppler width)
// Z1 = sqrt(X+Y) - sqrt(Y)
// Z2 = sqrt(X+Y) + sqrt(Y)
// vp = sqrt(2*kB*T/m), the most probable speed
// nuVC(v) is the  speed-dependent  velocity-changing collision frequency, given by:
//               nuVC(v) = nuVC - eta*[gamma(v)-i*delta(v)]
// where nuVC is the frequency of velocity-changing collisions when assuming no correlation between velocity and rotational-state changes
// Using quadratic dependences of the line width and shift on the speed v, i.e.:
///     gamma(v)-i*delta(v) = C0 + C2*{(v/vp)^2-1.5}
// C0Tilde = (1-eta)*(C0-3*C0/2) + nuVC
// C0 = gamma0 - i*delta0
// C2Tilde = (1 - eta)*C2
// C2 = (gamma2 - i*delta2)
//
// ============================================================================================================================================================================================================================ //
// References:
// [1] NH Ngo, D Lisak, H Tran, JM Hartmann. An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes
//     J Quant Spectrosc Radiat Transf 2013
// [2] Boone CD, Walker KA, Bernath PF. Speed-dependent Voigt profile for water vapor in infrared remote sensing applications.
//     J Quant Spectrosc Radiat Transf 2007;105:525-532.
//
// ============================================================================================================================================================================================================================ //

																				//************************************************************************************//
																				// PARTIALLY-CORRELATED QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE (BOONE MODEL) //
																				//************************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | partially Correlated quadratic Speed Dependent Hard Collision (Boone model)        |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |   Correlation parameter    |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//  In this case, the  complex  probability function w(z) is computed thanks to the CERNlib routine
// ============================================================================================================================================================================================================================ //

long simulation::doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double eta = (fabs(*(ptrParam+parameter_3))<=(double)1.0) ? fabs(*(ptrParam+parameter_3)):(double)1.0;
	if (!eta)
		return doQuadraticSpeedDependentNelkinGhatakBoone(area,y,ptrParam);

	std::vector<double> paramProfile(17,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double nuVC = fabs(*(ptrParam+velocityChanging));
	const double gamma0 = *(ptrParam+gammaHWHM);
	const double delta0 = *(ptrParam+parameter_0);
	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx ZERO(0.0,0.0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = ((double)1.0-eta)*(c0 - (double)1.5*c2)+nuVC;
	const dcplx c2Tilde = ((double)1.0-eta)*c2;
	dcplx invC2Tilde(0.0,0.0);
	dcplx yCplx(0.0,0.0); 
	dcplx y2Cplx(0.0,0.0);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	if (c2Tilde != ZERO)
	{
		invC2Tilde = fcnInverse(c2Tilde);
		yCplx = ((double)0.5*invC2Tilde)*gHWHM;
		y2Cplx = fcnCarre(yCplx);
	}
	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = c2Tilde.real();
	paramProfile[8] = c2Tilde.imag();
	paramProfile[9] = invC2Tilde.real();
	paramProfile[10] = invC2Tilde.imag();
	paramProfile[11] = nuVC;
	paramProfile[12] = gamma0;
	paramProfile[13] = -delta0;
	paramProfile[14] = gamma2;
	paramProfile[15] = -delta2;
	paramProfile[16] = eta;

	const double amp = (!area) ? (amplitude/partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //
//	"pcqSDNGP_Humlicek": partially-Correlated quadratic-Speed-Dependent Hard-Collision
//	Subroutine to Compute the complex normalized spectral shape of an
//	isolated line by the pCqSDHC model
//
//	Called Routines: 'complexErrFunc'	(complex error function~Complex Probability Function)
// ============================================================================================================================================================================================================================ //

double simulation::partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{
	static const dcplx iCplx((double)0.0,(double)1.0);
	static const dcplx ZEROCplx((double)0.0,(double)1.0);
	dcplx c2Tilde(*(ptrParamP+7),*(ptrParamP+8));
	dcplx invC2Tilde(*(ptrParamP+9),*(ptrParamP+10));
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx y2Cplx(*(ptrParamP+5),*(ptrParamP+6));
	dcplx c0(*(ptrParamP+12),*(ptrParamP+13));
	dcplx c2(*(ptrParamP+14),*(ptrParamP+15));
	dcplx Aterm(0.0,0.0);
	dcplx Bterm(0.0,0.0);
	dcplx qpcSDNGP(0.0,0.0);
	dcplx z1CPlx(0.0,0.0);
	dcplx z2CPlx(0.0,0.0);
	dcplx w1CPlx(0.0,0.0);
	dcplx w2CPlx(0.0,0.0);
	dcplx xCplx(0.0,0.0);
	double borneSup((double)3.0e-8*abs(y2Cplx));
	double center(0.0);
	double xz1(0.0);
	double xz2(0.0);
	double yz1(0.0);
	double yz2(0.0);
	double SZ1(0.0);
	double SZ2(0.0);
	double DSZ(0.0);
	double SZmax(0.0);
	double SZmin(0.0);
	long region(0);

	if (c2Tilde!=ZEROCplx)
		xCplx = ((iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*invC2Tilde);
	else
	{
		// when C2t=0
		region = 3;
	}
	// when abs(X) is much larger than abs(Y)
	if (abs(yCplx)<(double)1.0e-15*abs(xCplx))
		region = 1;
	// when abs(Y) is much larger than abs(X)
	else if (abs(xCplx)<borneSup)
		region = 2;
	switch(region)
	{
		case 0:{					
					// calculating Z1 and Z2
					z1CPlx = sqrt(xCplx+y2Cplx)-yCplx;
					z2CPlx = z1CPlx + (double)2.0*yCplx;
					// calculating the real and imaginary parts of Z1 and Z2
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					// check if Z1 and Z2 are close to each other
					SZ1 = fcnNorme(xz1,yz1);
					SZ2 = fcnNorme(xz2,yz2);
					DSZ = fabs(SZ1-SZ2);
					SZmax = max(SZ1,SZ2);
					SZmin = min(SZ1,SZ2);
					w1CPlx = complexErrFunc(dcplx(xz1,yz1));
					w2CPlx = complexErrFunc(dcplx(xz2,yz2));
					// calculating the A term of the profile
					Aterm = *(ptrParamP)*SQRTPI*(w1CPlx-w2CPlx);
					Bterm = (-(double)1.0 + SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z1CPlx))*w1CPlx
						                   -SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
			   }
			   break;
		case 1:{
					// when abs(X) is much larger than abs(Y)
					xz1 = -(sqrt(xCplx+y2Cplx)).imag();
					yz1 = (sqrt(xCplx+y2Cplx)).real();
					w1CPlx = complexErrFunc(dcplx(xz1,yz1));
					if (abs(sqrt(xCplx))<=(double)4000.0)
					{
						w2CPlx = complexErrFunc(dcplx(-(sqrt(xCplx)).imag(),sqrt(xCplx).real()));
						Aterm = ((double)2.0*SQRTPI*invC2Tilde)*(INVERSE_SQRTPI-sqrt(xCplx)*w1CPlx);
						Bterm = invC2Tilde*(-(double)1.0 + (double)2.0*SQRTPI*((double)1.0 - xCplx - (double)2.0*y2Cplx)*(INVERSE_SQRTPI-sqrt(xCplx)*w2CPlx)
							                        + (double)2.0*SQRTPI*sqrt(xCplx+y2Cplx)*w1CPlx);
					}
					else
					// when abs(X) is much larger than 1
					{
						Aterm = invC2Tilde*(fcnInverse(xCplx)-(double)1.5*fcnInverse(fcnCarre(xCplx)));
						Bterm = invC2Tilde*(-(double)1.0+((double)1.0 - xCplx - (double)2.0*y2Cplx)*(fcnInverse(xCplx)-(double)1.5*fcnInverse(fcnCarre(xCplx)))
							                                                                                          +(double)2.0*SQRTPI*sqrt(xCplx+y2Cplx)*w1CPlx);
					}
			   }
			   break;
		case 2:{
					// when abs(Y) is much larger than abs(X)
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = *(ptrParamP)*(iCplx*x+c0t);
					z2CPlx = sqrt(xCplx+y2Cplx) + yCplx;
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					w1CPlx = complexErrFunc(dcplx(xz1,yz1));
					w2CPlx = complexErrFunc(dcplx(xz2,yz2));
					Aterm = *(ptrParamP)*SQRTPI*(w1CPlx-w2CPlx);
					Bterm = (-(double)1.0 + SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z1CPlx))*w1CPlx
						                   -SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
			   }
			   break;
		case 3:{
					// when C2t=0
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = *(ptrParamP)*(iCplx*x+c0t);
					w1CPlx = complexErrFunc(dcplx(-z1CPlx.imag(),z1CPlx.real()));
					Aterm = *(ptrParamP)*SQRTPI*w1CPlx;
					Bterm = *(ptrParamP)*((SQRTPI*w1CPlx) + (double)0.5/z1CPlx - (double)0.75*Power<3>::of(z1CPlx));
			   }
			   break;

		default:break;
	}
	qpcSDNGP = INVERSE_PI*Aterm/((double)1.0-(*(ptrParamP+11)-*(ptrParamP+16)*(c0 - (double)1.5*c2))*Aterm + (*(ptrParamP+14))*c2*Bterm);

	return (qpcSDNGP.real());

}

// ============================================================================================================================================================================================================================ //

																				//**********************************************************************************************//
																				// AREA PARTIALLY-CORRELATED QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE (HUMLICEK ROUTINE) //
																				//**********************************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | partially Correlated quadratic Speed Dependent Hard Collision (Boone model)        |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |   Correlation parameter    |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double eta = (fabs(*(ptrParam+parameter_3))<=(double)1.0) ? fabs(*(ptrParam+parameter_3)):(double)1.0;
	if (!eta)
		return doAreaQuadraticSpeedDependentNelkinGhatakBoone(ptrParam);

	std::vector<double> paramProfile(17,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double nuVC = fabs(*(ptrParam+velocityChanging));
	const double gamma0 = *(ptrParam+gammaHWHM);
	const double delta0 = *(ptrParam+parameter_0);
	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx ZERO(0.0,0.0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = ((double)1.0-eta)*(c0 - (double)1.5*c2)+nuVC;
	const dcplx c2Tilde = ((double)1.0-eta)*c2;
	dcplx invC2Tilde(0.0,0.0);
	dcplx yCplx(0.0,0.0); 
	dcplx y2Cplx(0.0,0.0);
	double lineCenter((double)0.0);

	if (c2Tilde != ZERO)
	{
		invC2Tilde = fcnInverse(c2Tilde);
		yCplx = ((double)0.5*invC2Tilde)*gHWHM;
		y2Cplx = fcnCarre(yCplx);
	}
	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = c2Tilde.real();
	paramProfile[8] = c2Tilde.imag();
	paramProfile[9] = invC2Tilde.real();
	paramProfile[10] = invC2Tilde.imag();
	paramProfile[11] = nuVC;
	paramProfile[12] = gamma0;
	paramProfile[13] = -delta0;
	paramProfile[14] = gamma2;
	paramProfile[15] = -delta2;
	paramProfile[16] = eta;

	const double norm = fcnInverse(partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// In  Ref[2],  it  was  shown  that  the  quadratic-speed-dependent  Voigt  profile  can  be
// expressed as a combination of two "usual" Voigt functions.
// Within the  partially-Correlated Speed-Dependent Hard-Collision model of an isolated
// line is given by can be written as:
//      IpCqSDHC(u) = 1/PI*{A(u)/([1-nuVC-eta(C0-1.5*C2)]*A(u) + [eta*C2/vp^2)*B(u)])}
// The A(u) and B(u) terms can be expressed as combination of the complex probability function w(z):
// 1/PI*A(u) = c/(sqrt(PI)*u0*vp)*[w(iZ1)-w(iZ2)]
// B(u) = (vp^2/C2Tilde)*[-1+sqrt(PI)/(2*sqrt(Y)*{(1-Z1^2)*w(iZ1)-(1-Z2^2)*w(iZ2)})]
// with
// X = [i(u-u0)+C0Tilde]/C2Tilde
// Y = (u0*vp/(2*c*C2Tilde))^2 = (gammaD/(2*C2*sqrt(ln2)))^2 (gammaD is the Doppler width)
// Z1 = sqrt(X+Y) - sqrt(Y)
// Z2 = sqrt(X+Y) + sqrt(Y)
// vp = sqrt(2*kB*T/m), the most probable speed
// nuVC(v) is the  speed-dependent  velocity-changing collision frequency, given by:
//               nuVC(v) = nuVC - eta*[gamma(v)-i*delta(v)]
// where nuVC is the frequency of velocity-changing collisions when assuming no correlation between velocity and rotational-state changes
// Using quadratic dependences of the line width and shift on the speed v, i.e.:
///     gamma(v)-i*delta(v) = C0 + C2*{(v/vp)^2-1.5}
// C0Tilde = (1-eta)*(C0-3*C0/2) + nuVC
// C0 = gamma0 - i*delta0
// C2Tilde = (1 - eta)*C2
// C2 = (gamma2 - i*delta2)
//
// ============================================================================================================================================================================================================================ //
// References:
// [1] NH Ngo, D Lisak, H Tran, JM Hartmann. An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes
//     J Quant Spectrosc Radiat Transf 2013
// [2] Boone CD, Walker KA, Bernath PF. Speed-dependent Voigt profile for water vapor in infrared remote sensing applications.
//     J Quant Spectrosc Radiat Transf 2007;105:525-532.
//
// ============================================================================================================================================================================================================================ //

																				//**********************************************************************************************************//
																				// PARTIALLY-CORRELATED QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE (BOONE MODEL WITH HUMLICEK ROUTINE) //
																				//**********************************************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | partially Correlated quadratic Speed Dependent Hard Collision                      |                                 |                            |                            |                            |
// |                       (Boone model with Humlicek routine)                          |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |   Correlation parameter    |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //
// Based on a slightly improved version of the CPF subroutine [Humlicek,J Quant Spectrosc Radiat Transf 21 309 (1979)]
// for the calculation of the complex probability function
// ============================================================================================================================================================================================================================ //

long simulation::doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const double eta = (fabs(*(ptrParam+parameter_3))<=(double)1.0) ? fabs(*(ptrParam+parameter_3)):(double)1.0;
	if (!eta)
		return doQuadraticSpeedDependentNelkinGhatakHumlicekBoone(area,y,ptrParam);

	std::vector<double> paramProfile(17,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double nuVC = fabs(*(ptrParam+velocityChanging));
	const double gamma0 = *(ptrParam+gammaHWHM);
	const double delta0 = *(ptrParam+parameter_0);
	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx ZERO(0.0,0.0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = ((double)1.0-eta)*(c0 - (double)1.5*c2)+nuVC;
	const dcplx c2Tilde = ((double)1.0-eta)*c2;
	dcplx invC2Tilde(0.0,0.0);
	dcplx yCplx(0.0,0.0); 
	dcplx y2Cplx(0.0,0.0);
	double lineCenter((double)0.0);
	double *ptrY = &y[0];

	if (c2Tilde != ZERO)
	{
		invC2Tilde = fcnInverse(c2Tilde);
		yCplx = ((double)0.5*invC2Tilde)*gHWHM;
		y2Cplx = fcnCarre(yCplx);
	}
	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = c2Tilde.real();
	paramProfile[8] = c2Tilde.imag();
	paramProfile[9] = invC2Tilde.real();
	paramProfile[10] = invC2Tilde.imag();
	paramProfile[11] = nuVC;
	paramProfile[12] = gamma0;
	paramProfile[13] = -delta0;
	paramProfile[14] = gamma2;
	paramProfile[15] = -delta2;
	paramProfile[16] = eta;

	const double amp = (!area) ? (amplitude/partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center);
		(*ptrY++)+=amp*partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //
//	"pcqSDNGP_Humlicek": partially-Correlated quadratic-Speed-Dependent Hard-Collision
//	Subroutine to Compute the complex normalized spectral shape of an
//	isolated line by the pCqSDHC model
//
//	Called Routines: 'CPF'	(Complex Probability Function)
//	---------------  'CPF3'	(Complex Probability Function for the region 3)
// ============================================================================================================================================================================================================================ //

double simulation::partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	static const dcplx iCplx((double)0.0,(double)1.0);
	static const dcplx ZEROCplx((double)0.0,(double)1.0);
	dcplx c2Tilde(*(ptrParamP+7),*(ptrParamP+8));
	dcplx invC2Tilde(*(ptrParamP+9),*(ptrParamP+10));
	dcplx yCplx(*(ptrParamP+3),*(ptrParamP+4));
	dcplx y2Cplx(*(ptrParamP+5),*(ptrParamP+6));
	dcplx c0(*(ptrParamP+12),*(ptrParamP+13));
	dcplx c2(*(ptrParamP+14),*(ptrParamP+15));
	dcplx Aterm(0.0,0.0);
	dcplx Bterm(0.0,0.0);
	dcplx qpcSDNGP(0.0,0.0);
	dcplx z1CPlx(0.0,0.0);
	dcplx z2CPlx(0.0,0.0);
	dcplx w1CPlx(0.0,0.0);
	dcplx w2CPlx(0.0,0.0);
	dcplx xCplx(0.0,0.0);
	double borneSup((double)3.0e-8*abs(y2Cplx));
	double center(0.0);
	double xz1(0.0);
	double xz2(0.0);
	double yz1(0.0);
	double yz2(0.0);
	double SZ1(0.0);
	double SZ2(0.0);
	double DSZ(0.0);
	double SZmax(0.0);
	double SZmin(0.0);
	long region(0);

	if (c2Tilde!=ZEROCplx)
		xCplx = ((iCplx*x+dcplx(*(ptrParamP+1),*(ptrParamP+2)))*invC2Tilde);
	else
	{
		// when C2t=0
		region = 3;
	}
	// when abs(X) is much larger than abs(Y)
	if (abs(yCplx)<(double)1.0e-15*abs(xCplx))
		region = 1;
	// when abs(Y) is much larger than abs(X)
	else if (abs(xCplx)<borneSup)
		region = 2;
	switch(region)
	{
		case 0:{
					invC2Tilde = fcnInverse(c2Tilde);
					// calculating Z1 and Z2
					z1CPlx = sqrt(xCplx+y2Cplx)-yCplx;
					z2CPlx = z1CPlx + (double)2.0*yCplx;
					// calculating the real and imaginary parts of Z1 and Z2
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					// check if Z1 and Z2 are close to each other
					SZ1 = fcnNorme(xz1,yz1);
					SZ2 = fcnNorme(xz2,yz2);
					DSZ = fabs(SZ1-SZ2);
					SZmax = max(SZ1,SZ2);
					SZmin = min(SZ1,SZ2);
					// when Z1 and Z2 are close to each other, ensure that they are in
					// the same interval of CPF
					if (DSZ<=(double)1.0 && SZmax>(double)8.0 &&  SZmin<=(double)8.0)
					{
						w1CPlx = CPF3(dcplx(xz1,yz1));
						w2CPlx = CPF3(dcplx(xz2,yz2));
					}
					else
					{
						w1CPlx = CPF(dcplx(xz1,yz1));
						w2CPlx = CPF(dcplx(xz2,yz2));
					}
					// calculating the A term of the profile
					Aterm = *(ptrParamP)*SQRTPI*(w1CPlx-w2CPlx);
					Bterm = (-(double)1.0 + SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z1CPlx))*w1CPlx
						                   -SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
			   }
			   break;
		case 1:{
					invC2Tilde = fcnInverse(c2Tilde);
					xz1 = -(sqrt(xCplx+y2Cplx)).imag();
					yz1 = (sqrt(xCplx+y2Cplx)).real();
					w1CPlx = CPF(dcplx(xz1,yz1));
					// when abs(X) is much larger than abs(Y)
					if (abs(sqrt(xCplx))<=(double)4000.0)
					{
						w2CPlx = CPF(dcplx(-(sqrt(xCplx)).imag(),sqrt(xCplx).real()));
						Aterm = ((double)2.0*SQRTPI*invC2Tilde)*(INVERSE_SQRTPI-sqrt(xCplx)*w1CPlx);
						Bterm = invC2Tilde*(-(double)1.0 + (double)2.0*SQRTPI*((double)1.0 - xCplx - (double)2.0*y2Cplx)*(INVERSE_SQRTPI-sqrt(xCplx)*w2CPlx)
							                             + (double)2.0*SQRTPI*sqrt(xCplx+y2Cplx)*w1CPlx);
					}
					else
					// when abs(X) is much larger than 1
					{
						Aterm = invC2Tilde*(fcnInverse(xCplx)-(double)1.5*fcnInverse(fcnCarre(xCplx)));
						Bterm = invC2Tilde*(-(double)1.0+((double)1.0 - xCplx - (double)2.0*y2Cplx)*(fcnInverse(xCplx)-(double)1.5*fcnInverse(fcnCarre(xCplx)))
							                                                                                          +(double)2.0*SQRTPI*sqrt(xCplx+y2Cplx)*w1CPlx);
					}
			   }
			   break;
		case 2:{
					invC2Tilde = fcnInverse(c2Tilde);
					// when abs(Y) is much larger than abs(X)
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = *(ptrParamP)*(iCplx*x+c0t);
					z2CPlx = sqrt(xCplx+y2Cplx) + yCplx;
					xz1 = -z1CPlx.imag();
					yz1 = z1CPlx.real();
					xz2 = -z2CPlx.imag();
					yz2 = z2CPlx.real();
					w1CPlx = CPF(dcplx(xz1,yz1));
					w2CPlx = CPF(dcplx(xz2,yz2));
					Aterm = *(ptrParamP)*SQRTPI*(w1CPlx-w2CPlx);
					Bterm = (-(double)1.0 + SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z1CPlx))*w1CPlx
						                   -SQRTPI/((double)2.0*yCplx)*((double)1.0-fcnCarre(z2CPlx))*w2CPlx)*invC2Tilde;
			   }
			   break;
		case 3:{
					// when C2t=0
					dcplx c0t(*(ptrParamP+1),*(ptrParamP+2));
					z1CPlx = *(ptrParamP)*(iCplx*x+c0t);
					w1CPlx = CPF(dcplx(-z1CPlx.imag(),z1CPlx.real()));
					Aterm = *(ptrParamP)*SQRTPI*w1CPlx;
					if (abs(sqrt(z1CPlx))<(double)4000.0)
					{
						Bterm = *(ptrParamP)*SQRTPI*((((double)1.0-fcnCarre(z1CPlx))*w1CPlx) + *(ptrParamP)*z1CPlx);
					}
					else
					{
						Bterm = *(ptrParamP)*((SQRTPI*w1CPlx) + (double)0.5/z1CPlx - (double)0.75*Power<3>::of(z1CPlx));
					}
			   }
			   break;
		default:break;
	}
	qpcSDNGP = INVERSE_PI*Aterm/((double)1.0-(*(ptrParamP+11)-*(ptrParamP+16)*(c0 - (double)1.5*c2))*Aterm + (*(ptrParamP+14))*c2*Bterm);

	return (qpcSDNGP.real());

}

// ============================================================================================================================================================================================================================ //

																				//*****************************************************************************************//
																				// PARTIALLY-CORRELATED QUADRATIC SPEED-DEPENDENT NELKIN-GHATAK PROFILE (HUMLICEK ROUTINE) //
																				//*****************************************************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | partially Correlated quadratic Speed Dependent Hard Collision (Humlicek)	        |	         Shift #0	          |	        Gamma #2	       |  		 Shift #2           |   Correlation parameter    |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone(std::vector<double>::const_iterator ptrParam)
{
	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const double eta = (fabs(*(ptrParam+parameter_3))<=(double)1.0) ? fabs(*(ptrParam+parameter_3)):(double)1.0;
	if (!eta)
		return doAreaQuadraticSpeedDependentNelkinGhatakHumlicekBoone(ptrParam);

	std::vector<double> paramProfile(17,(double)0.0);
	std::vector<double>::const_iterator ptrParamProfile = paramProfile.begin();
	const double nuVC = fabs(*(ptrParam+velocityChanging));
	const double gamma0 = *(ptrParam+gammaHWHM);
	const double delta0 = *(ptrParam+parameter_0);
	const double gamma2 = fabs(*(ptrParam+parameter_1));
	const double delta2 = *(ptrParam+parameter_2);
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const dcplx ZERO(0.0,0.0);
	const dcplx c0(gamma0,-delta0);
	const dcplx c2(gamma2,-delta2);
	const dcplx c0Tilde = ((double)1.0-eta)*(c0 - (double)1.5*c2)+nuVC;
	const dcplx c2Tilde = ((double)1.0-eta)*c2;
	dcplx invC2Tilde(0.0,0.0);
	dcplx yCplx(0.0,0.0); 
	dcplx y2Cplx(0.0,0.0);
	double lineCenter((double)0.0);

	if (c2Tilde != ZERO)
	{
		invC2Tilde = fcnInverse(c2Tilde);
		yCplx = ((double)0.5*invC2Tilde)*gHWHM;
		y2Cplx = fcnCarre(yCplx);
	}
	paramProfile[0] = invGWidth;
	paramProfile[1] = c0Tilde.real();
	paramProfile[2] = c0Tilde.imag();
	paramProfile[3] = yCplx.real();
	paramProfile[4] = yCplx.imag();
	paramProfile[5] = y2Cplx.real();
	paramProfile[6] = y2Cplx.imag();
	paramProfile[7] = c2Tilde.real();
	paramProfile[8] = c2Tilde.imag();
	paramProfile[9] = invC2Tilde.real();
	paramProfile[10] = invC2Tilde.imag();
	paramProfile[11] = nuVC;
	paramProfile[12] = gamma0;
	paramProfile[13] = -delta0;
	paramProfile[14] = gamma2;
	paramProfile[15] = -delta2;
	paramProfile[16] = eta;

	const double norm = fcnInverse(partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);

}

// ============================================================================================================================================================================================================================ //
//
// The approximate speed-dependent Galatry Profile (asdGP) is expressed as:
//
// IasdGP(x,y,z) = 1/PI*{$$\int_(0^+INF)exp[-(1/4z^2)*(2*z*t-3+4*exp(-zt)-exp(-2*z*t)]
//            *4/sqrt(PI)*[\int_(0^+INF)sinc[(1-exp(-z*t))*u/z]*exp{-[ix(u)+y(u)]t}exp(-u^2)u^2du]dt$$}
//
// Assuming a quadratic  dependence  for  each  of  the  speed-dependent
// collisional functions then, in dimensionless variable:
//
//                y(u) = y*[1 + aw*(u^2 - 3/2)],
//				  x(u) = x - s(u) with s(u) = s*[1 + as*(u^2 - 3/2)].
//
// ============================================================================================================================================================================================================================ //
// Reference:
// Collisional Effects on Molecular Spectra
// Laboratory Experiments and Models, Consequences for Applications
// Author(s): Jean-Michel Hartmann, Christian Boulet and Daniel Robert
// ============================================================================================================================================================================================================================ //

																				//*******************************************//
																				// QUADRATIC SPEED-DEPENDENT GALATRY PROFILE //
																				//*******************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | approximate quadratic Speed Dependent Galatry				                        |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

long simulation::doApproximateQuadraticSpeedDependentGalatry(long &area, double *y, std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return OK;
	if(amplitude<(double)0.0)
		return ERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return ERREUR;

	const long SDGProfile(7); 
	const long nbParam = SDGProfile + 3*NB_PTS_GAUSS_LEGENDRE_SDG_INT;
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double zs = fabs((*(ptrParam+velocityChanging))*invGWidth);
	double *ptrY = &y[0];
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = FOURbyCUBERTPI*invGWidth;
	paramProfile[1] = fabs((*(ptrParam+gammaHWHM))*invGWidth);
	paramProfile[2] = (*(ptrParam+parameter_0))*invGWidth;
	paramProfile[3] = (*(ptrParam+parameter_1));
	paramProfile[4] = (*(ptrParam+parameter_2));
	paramProfile[5] = zs;
	paramProfile[6] = (double)0.25/fcnCarre(zs);
	initQSDGP(SDGProfile,ptrParamProfile); 

	const double amp = (!area) ? (amplitude/approximateQuadraticSpeedDependentGalatryProfile(ptrParamProfile,lineCenter)):amplitude;

	for(std::vector<double>::const_iterator ptrFx(fX.begin()),ptrFxEnd(fX.end());ptrFx!=ptrFxEnd;++ptrFx)
	{
		lineCenter = ((*ptrFx)-center)*invGWidth;
		(*ptrY++)+=amp*approximateQuadraticSpeedDependentGalatryProfile(ptrParamProfile,lineCenter);
	}

	return OK;
}

// ============================================================================================================================================================================================================================ //

double simulation::approximateQuadraticSpeedDependentGalatryProfile(std::vector<double>::const_iterator ptrParamP,double &x)
{

	const long shiftSDG = 7;
	double argTrig(0.0);
	double dummy(0.0);
	double dummyZero(0.0);
	std::vector<double>::const_iterator integraleSDGPExp(ptrParamP+shiftSDG);
	std::vector<double>::const_iterator integraleSDGPCos(integraleSDGPExp+NB_PTS_GAUSS_LEGENDRE_SDG_INT);
	std::vector<double>::const_iterator integraleSDGPSin(integraleSDGPCos+NB_PTS_GAUSS_LEGENDRE_SDG_INT);

	for(long j=0;j<NB_PTS_GAUSS_LEGENDRE_SDG_INT;++j)
	{
		argTrig = x*xqSDG4096[j];
		dummy=(*integraleSDGPExp++)*(cos(argTrig)*(*integraleSDGPCos++) + sin(argTrig)*(*integraleSDGPSin++));
		if (isnan(dummy)|| isinf(dummy))
			break;	
		dummyZero +=dummy;	
	}

	return (*(ptrParamP)*dummyZero);

}

// ============================================================================================================================================================================================================================ //

																				//************************************************//
																				// AREA QUADRATIC SPEED-DEPENDENT GALATRY PROFILE //
																				//************************************************//

// ============================================================================================================================================================================================================================ //
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
//                                PROFILE\PARAMETER		   				     	        |	            #0	              |	          #1	   	       | 		     #2	            | 	      #3    	         |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// | approximate quadratic Speed Dependent Galatry				                        |	           Shift 	     	  |		      aw     	       |		    as              |	     		             |
// -------------------------------------------------------------------------------------|---------------------------------|----------------------------|----------------------------|----------------------------|
// ============================================================================================================================================================================================================================ //

double simulation::doAreaApproximateQuadraticSpeedDependentGalatry(std::vector<double>::const_iterator ptrParam)
{

	const double amplitude = *(ptrParam+Amplitude);
	if(!amplitude)
		return dblOK;
	if(amplitude<(double)0.0)
		return dblERREUR;

	const double gHWHM = fabs(*(ptrParam+dopplerHWHM));
	if (!gHWHM)
		return dblERREUR;

	const long SDGProfile(7); 
	const long nbParam = SDGProfile + 3*NB_PTS_GAUSS_LEGENDRE_SDG_INT;
	const double invGWidth = fcnInverse(gHWHM);
	const double center = *(ptrParam+Wavelength);
	const double zs = fabs(*(ptrParam+velocityChanging))*invGWidth;
	double lineCenter((double)0.0);
	std::vector<double> paramProfile(nbParam,(double)0.0);
	std::vector<double>::iterator ptrParamProfile = paramProfile.begin();

	paramProfile[0] = FOURbyCUBERTPI*invGWidth;
	paramProfile[1] = fabs((*(ptrParam+gammaHWHM))*invGWidth);
	paramProfile[2] = (*(ptrParam+parameter_0))*invGWidth;
	paramProfile[3] = (*(ptrParam+parameter_1));
	paramProfile[4] = (*(ptrParam+parameter_2));
	paramProfile[5] = zs;
	paramProfile[6] = (double)0.25/fcnCarre(zs);
	initQSDGP(SDGProfile,ptrParamProfile);
	
	const double norm = fcnInverse(approximateQuadraticSpeedDependentGalatryProfile(ptrParamProfile,lineCenter));

	return (amplitude*norm);
}

// ============================================================================================================================================================================================================================ //

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																		//************************************//
																		// Pour le calcul de la ligne de base //
																		//************************************//

// ============================================================================================================================================================================================================================ //

// ********************************************* //
// Calcul des abscisses [-1;1] pour le polynôme  //
// MàJ si le nombre des données change           //
// ********************************************* //

inline void simulation::initAbscisse(std::vector<double> &abscisse)
{
	std::vector<double>::const_iterator ptrFX(fX.begin());
	std::vector<double>::iterator ptrAbs(abscisse.begin());
	double nbData  = (double)fX.size();
	double index(0.0);
	double a = (*ptrFX);
	double b = *(fX.end()-1);
	double c(0.0);
	double d(0.0);

	if (a>b)
	{
		double temp(0.0);
		temp = a;
		a = b;
		b = temp;
	}
	c = (double)2.0/(b-a);
	d = (double)0.5*(b+a);
	for (; ptrFX != fX.end(); ++ptrFX,++ptrAbs)
	{
		(*ptrAbs) = c*((*ptrFX)-d);
	}
}

// ============================================================================================================================================================================================================================ //

void simulation::FCN_Baseline(double *y, const std::vector<double> &fPar)
{
	const long nbCoeffsPolynome=(long)fPar[1];
	const long nbData  = (long)fX.size();
	const double begin(fX[0]);
	const double end(fX[nbData-1]);
	static double oldBegin(0.0);
	static double oldEnd(0.0);
	static double oldNbData(0.0);
	static std::vector<double> abscisseBaseline(nbData,0.0);
	bool firstCall = ((oldBegin != begin)|| (oldEnd != end) || (oldNbData != nbData));

	if (firstCall && nbCoeffsPolynome>1)
	{
		abscisseBaseline.clear();
		if (nbData!=abscisseBaseline.size())
		{
			oldNbData = nbData;
			abscisseBaseline.resize(nbData,0.0);
		}
		oldBegin = begin;
		oldEnd = end;
		initAbscisse(abscisseBaseline);
	}
	switch (nbCoeffsPolynome)
	{
		case 1 :{
					BaselineConstant(y,fPar);
				}
				break;
		case 2 :{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomial1(y,fPar,ptrAbs);
				}
				break;
		case 3 :{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomial2(y,fPar,ptrAbs);
				}
				break;
		case 4 :{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomial3(y,fPar,ptrAbs);
				}
				break;
		case 5 :{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomial4(y,fPar,ptrAbs);
				}
				break;
		case 6 :{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomial5(y,fPar,ptrAbs);
				}
				break;
		default:{
					std::vector<double>::const_iterator ptrAbs(abscisseBaseline.begin());
					BaselinePolynomialn(y,fPar,ptrAbs);
				}
				break;
	}
}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 0 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselineConstant(double *y, const std::vector<double> &coeffBaseline)
{

	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;

	a0 = *argumentX;
	for (long i=0;i<nbData;++i)
	{
		(*y)+=a0;
		++y;
	}

}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 1 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomial1(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	double a1(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);

	a0 = *argumentX;
	a1 = *(++argumentX);
	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		(*y)+=a1*(*ptrAbs)+a0;
		++y;
	}

}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 2 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomial2(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	double a1(0.0);
	double a2(0.0);
	double x(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);

	a0 = *argumentX;
	a1 = *(++argumentX);
	a2 = *(++argumentX);
	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		x = (*ptrAbs);
		(*y)+= Power<2>::of(x)*a2 + x*a1 +a0;
		++y;
	}
}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 3 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomial3(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	double a1(0.0);
	double a2(0.0);
	double a3(0.0);
	double x(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);

	a0 = *argumentX;
	a1 = *(++argumentX);
	a2 = *(++argumentX);
	a3 = *(++argumentX);
	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		x = (*ptrAbs);
		(*y)+= Power<3>::of(x)*a3 + Power<2>::of(x)*a2 + x*a1 +a0;
		++y;
	}
}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 4 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomial4(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	double a1(0.0);
	double a2(0.0);
	double a3(0.0);
	double a4(0.0);
	double x(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);

	a0 = *argumentX;
	a1 = *(++argumentX);
	a2 = *(++argumentX);
	a3 = *(++argumentX);
	a4 = *(++argumentX);
	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		x = (*ptrAbs);
		(*y)+= Power<4>::of(x)*a4 + Power<3>::of(x)*a3 + Power<2>::of(x)*a2 + x*a1 +a0;
		++y;
	}
}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree 5 polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomial5(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	long nbData  = (long)fX.size();
	long nbCoeffPolynome = (long)fPar[1];
	double a0(0.0);
	double a1(0.0);
	double a2(0.0);
	double a3(0.0);
	double a4(0.0);
	double a5(0.0);
	double x(0.0);
	std::vector<double>::const_iterator argumentX = coeffBaseline.end()- nbCoeffPolynome;
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);

	a0 = *argumentX;
	a1 = *(++argumentX);
	a2 = *(++argumentX);
	a3 = *(++argumentX);
	a4 = *(++argumentX);
	a5 = *(++argumentX);
	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		x = (*ptrAbs);
		(*y)+= Power<5>::of(x)*a5 +Power<4>::of(x)*a4 + Power<3>::of(x)*a3 + Power<2>::of(x)*a2 + x*a1 +a0;
		++y;
	}
}

// ============================================================================================================================================================================================================================ //

																		//*********************//
																		// degree n polynomial //
																		//*********************//

// ============================================================================================================================================================================================================================ //

void simulation::BaselinePolynomialn(double *y, const std::vector<double> &coeffBaseline, std::vector<double>::const_iterator ptrAbs)
{
	const long nbData  = (long)fX.size();
	const long nbCoeffPolynome = (long)fPar[1];
	const long degree = nbCoeffPolynome -1;
	std::vector<double>::const_iterator first = coeffBaseline.end()- nbCoeffPolynome;
    std::vector<double>::const_iterator last  = coeffBaseline.end();
	std::vector<double> polynomeCoeff(first,last);
	std::vector<double>::const_iterator ptrAbsEnd(ptrAbs+nbData);
	double x(0.0);
	double an = *(--last);
	double sum(0.0);

	for (; ptrAbs != ptrAbsEnd; ++ptrAbs)
	{
		x = (*ptrAbs);
		sum = an;
		for(long j = degree - 1; j >= 0; --j)
		{
			sum *= x;
			sum += polynomeCoeff[j];
		}
		(*y)+=sum;
		++y;
	}
}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																				// ********************************************** //
																				/*         Fonctions auxiliaires                  */
																				// ********************************************** //

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //


// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

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

dcplx simulation::complexErrFunc(const dcplx &z)
{
  // ======================================================================================== //
  /* Return CERNlib complex error function                                                    */
  /*                                                                                          */
  /* This code is translated from the fortran version in the CERN mathlib.                    */
  /* (see http://aliceinfo.cern.ch/alicvs/viewvc/MINICERN/mathlib/gen/c/cwerf64.F?view=co)    */
  // ======================================================================================== //

    static std::vector<dcplx> r(38);
								// { HF, C1, C2,    C3, C4,                                P}
	static double lookup_table[] = {0.5,7.4,8.3,0.3125,1.6,46768052394588893.382517914646921};
    dcplx zh,s,t,v;
    static const dcplx zero(0.0,0.0);
    double x(z.real()),y(z.imag()), xAbs(fabs(x)), yAbs(fabs(y));
    int N;

    if((yAbs < lookup_table[1]) && (xAbs <  lookup_table[2]))
	{
        zh = dcplx(yAbs+ lookup_table[4],xAbs);
        r[37] = zero;
        N = 36;
        while(N > 0)
		{
            t = zh + conj(r[N+1])*(double)N;
            r[N--] = (t*lookup_table[0])/norm(t);
        }
        double xl = lookup_table[5];
        s = zero;
        N = 33;
        while(N > 0)
		{
            xl = lookup_table[3]*xl;
            s = r[N--]*(s+xl);
        }
        v = s*TWObySQRTPI;
    }
    else
	{
        zh = dcplx(yAbs,xAbs);
        r[1] = zero;
        N = 9;
        while(N > 0)
		{
            t = zh + conj(r[1])*(double)N;
            r[1] = (t*lookup_table[0])/norm(t);
            N--;
        }
        v = r[1]*TWObySQRTPI;
    }
    if (yAbs == (double)0.0)
		v = dcplx(exp(-(xAbs*xAbs)),v.imag());
    if (y < (double)0.0)
	{
        dcplx tmp(yAbs*yAbs - xAbs*xAbs,-(double)2.0*xAbs*yAbs);
        v = (double)2.0*exp(tmp) - v;
        if (x > (double)0.0)
			v = conj(v);
    }
    else
	{
        if(x < (double)0.0)
			v = conj(v);
    }
	return v;
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

dcplx simulation::CPF(const dcplx &z)
{

	dcplx  result(0.0,0.0);
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
	if(abs(z)>(double)8.0)
		region = 3;
	else if (y>(double)0.85||fabs(x)<((double)18.1*y+(double)1.65))
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

					for (long i=0;i<6;++i)
					{
						R = x - T[i];
						R2 = R*R;
						D=(double)1.0/(R2+Y2);
						D1 = Y1*D;
						D2 = R*D;
						R = x + T[i];
						R2 = R*R;
						D = (double)1.0/(R2+Y2);
						D3 = Y1*D;
						D4 = R*D;
						wr+=(U[i]*(D1+D3))-(S[i]*(D2-D4));
						wi+=(U[i]*(D2+D4))+(S[i]*(D1-D3));
					}
					result = dcplx(wr,wi);
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
					for (long i=0;i<6;++i)
					{
						R = x - T[i];
						R2 = R*R;
						D=(double)1.0/(R2+Y2);
						D1 = Y1*D;
						D2 = R*D;
						wr+=y*(U[i]*(R*D2-(double)1.5*D1)+S[i]*Y3*D2)/(R2+(double)2.25);
						R = x + T[i];
						R2 = R*R;
						D=(double)1.0/(R2+Y2);
						D3 = Y1*D;
						D4 = R*D;
						wr+= y*(U[i]*(R*D4-(double)1.5*D3)-S[i]*Y3*D4)/(R2+(double)2.25);
						wi+= U[i]*(D2+D4)+S[i]*(D1-D3);
					}
					result = dcplx(wr,wi);
				}
				break;
		case 3: {
					const dcplx zone(1.0,0.0);
					const dcplx zi(0.0,1.0);
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
					dcplx zm1(0.0,0.0);
					dcplx zm2(0.0,0.0);
					dcplx zterm(0.0,0.0);
					dcplx zsum(0.0,0.0);

					zm1=zone/z;
					zm2=zm1*zm1;
					zsum=zone;
					zterm=zone;
					for (long i=0;i<15;++i)
					{
						zterm*=zm2*tt[i];
						zsum+=zterm;
					}
					zsum*=zi*zm1*INVERSE_SQRTPI;
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

dcplx simulation::CPF3(const dcplx &z)
{
	dcplx zone(1.0,0.0);
	dcplx zi(0.0,1.0);
	dcplx zm1(0.0,0.0);
	dcplx zm2(0.0,0.0);
	dcplx zterm(0.0,0.0);
	dcplx zsum(0.0,0.0);
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

	zm1=zone/z;
	zm2=zm1*zm1;
	zsum=zone;
	zterm=zone;
	for (int i=0;i<15;++i)
	{
		zterm*=zm2*tt[i];
		zsum+=zterm;
	}
	zsum*=zi*zm1*pipwoeronehalf;

	return zsum;
}

// ============================================================================================================================================================== //

							// ***************************************************************************************************** //
							/*    Calcul des coefficients de l'expansion asymptotique pour la fonction GALATRY                       */
							// ***************************************************************************************************** //

// ============================================================================================================================================================== //

void simulation::coeffAsymtoticExpansion(std::vector<double>::iterator ptrCoefficientGalatry)
{

		const double z = *(ptrCoefficientGalatry+2);	
		double indice((double)4.0);
		ptrCoefficientGalatry[6] = (double)1.;
		ptrCoefficientGalatry[9] = z/(double)12.0;
		for(long j=10;j<15;++j)
		{
			ptrCoefficientGalatry[j]=-ptrCoefficientGalatry[j-1]*z/indice;
			indice+=(double)1.0;
		}
		ptrCoefficientGalatry[12]+=(double)0.5*fcnCarre(ptrCoefficientGalatry[9]);
		ptrCoefficientGalatry[13]+=ptrCoefficientGalatry[9]*ptrCoefficientGalatry[10];
		ptrCoefficientGalatry[14]+=ptrCoefficientGalatry[9]*ptrCoefficientGalatry[11]+(double)0.5*fcnCarre(ptrCoefficientGalatry[10]);
}

// ======================================================================================== //
// Initialisation des valeurs pour le calcul de l'intégrale dans le cas 
// "quadratic speed-dependent"
// ======================================================================================== //

void simulation::initQSDVP(const long &shiftProfile,std::vector<double>::iterator ptrParamP)
{
	std::vector<double>::iterator ptrBs(ptrParamP+shiftProfile);
	std::vector<double>::iterator ptrBw(ptrBs+NB_PTS_GAUSS_LEGENDRE);
	std::vector<double>::iterator ptrIntFcn(ptrBw+NB_PTS_GAUSS_LEGENDRE);
	double xi2(0.0);
	double xGL(0.0);

	for(long index =0;index<NB_PTS_GAUSS_LEGENDRE;++index)
	{
		xGL = PIby2*x1024[index];
		xi2 = tan(xGL)*tan(xGL);
		(*ptrBw++) = (*(ptrParamP+1))*((double)1.0 + *(ptrParamP+3)*(xi2-(double)1.5));
		(*ptrBs++) = (*(ptrParamP+2))*((double)1.0 + *(ptrParamP+4)*(xi2-(double)1.5));
		(*ptrIntFcn++) = w1024[index]*tan(xGL)*exp(-xi2)*(1.0+xi2);
	}

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void simulation::initQSDNGP(const long &shiftProfile,std::vector<double>::iterator ptrParamP)
{
	std::vector<double>::iterator ptrBs(ptrParamP+shiftProfile);
	std::vector<double>::iterator ptrBw(ptrBs+NB_PTS_GAUSS_LEGENDRE);
	std::vector<double>::iterator ptrIntFcn(ptrBw+NB_PTS_GAUSS_LEGENDRE);
	double xi2(0.0);
	double xGL(0.0);

	for(long index =0;index<NB_PTS_GAUSS_LEGENDRE;++index)
	{
		xGL = PIby2*x1024[index];
		xi2 = tan(xGL)*tan(xGL);
		(*ptrBw++) = (*(ptrParamP+1))*((double)1.0 + *(ptrParamP+3)*(xi2-(double)1.5)) + (*(ptrParamP+5));
		(*ptrBs++) = (*(ptrParamP+2))*((double)1.0 + *(ptrParamP+4)*(xi2-(double)1.5));
		(*ptrIntFcn++) = w1024[index]*tan(xGL)*exp(-xi2)*(1.0+xi2);
	}

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void simulation::initQSDGP(const long &shiftProfile,std::vector<double>::iterator ptrParamP)
{

//xqSDG4096
//dxqSDG4096
	
	std::vector<double>::iterator exponentielle(ptrParamP+shiftProfile);
	std::vector<double>::iterator integraleSDGPCos(exponentielle+NB_PTS_GAUSS_LEGENDRE_SDG_INT);
	std::vector<double>::iterator integraleSDGPSin(integraleSDGPCos+NB_PTS_GAUSS_LEGENDRE_SDG_INT);
	const double g(*(ptrParamP+1));
	const double s(*(ptrParamP+2));
	const double aw(*(ptrParamP+3));
	const double as(*(ptrParamP+4));
	const double zs(*(ptrParamP+5));
	const double invZs2(*(ptrParamP+6));
	const double invZs(fcnInverse(zs));
	double invU(0.0);
	double u2(0.0);
	double xGL(0.0);
	double argExp(0.0);
	double argExponentielle(0.0);
	double fctQD_Bw(0.0);
	double fctQD_Bs(0.0);
	double argTrig(0.0);
	double dummy(0.0);
	double argSinc(0.0);
	double dummyCos(0.0);
	double dummySin(0.0);


	for(long i=0;i<NB_PTS_GAUSS_LEGENDRE_SDG_INT;++i)
	{
		argExp = zs*xqSDG4096[i];
		argExponentielle = invZs2*((double)2.0*argExp-(double)3.0+(double)4.0*exp(-argExp)-exp(-(double)2.0*argExp));
		(*exponentielle++) = dxSDG4096[i]*exp(-argExponentielle);
		argSinc = ((double)1.0-exp(-argExp))*invZs;
		for(long j=0;j<NB_PTS_GAUSS_LEGENDRE_SDG_INT;++j)
		{
			u2 = fcnCarre(xqSDG4096[j]);
			fctQD_Bw = g*(aw*(u2-(double)1.5)+(double)1.);
			fctQD_Bs = s*(as*(u2-(double)1.5)+(double)1.);
			argTrig = fctQD_Bs*xqSDG4096[i];
			argExp = fctQD_Bw*xqSDG4096[i];
			dummy = dxSDG4096[j]*u2*exp(-u2)*fcnSinc(xqSDG4096[j]*argSinc)*exp(-argExp);
			if (isnan(dummy)|| isinf(dummy))
				break;	
			dummyCos +=dummy*cos(argTrig);
			dummySin +=dummy*sin(argTrig);
		}
		(*integraleSDGPCos++) = dummyCos;
		(*integraleSDGPSin++) = dummySin;
		dummyCos = (double)0.0;
		dummySin = (double)0.0;
	}

}

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

	}  // namespace Minuit2
} // namespace ROOT
