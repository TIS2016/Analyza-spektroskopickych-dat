#ifndef MN_SIMULATION_H_
#define MN_SIMULATION_H_

#include "mathFonctions.h"
#include "constantesProfils.h"
#include <complex>
#include <vector>

/************************************************************************************/
/* Structure Labview pour récupérer la fonction de covariance calculée par Minuit2  */
/************************************************************************************/

typedef struct {
	long dimSizes[2];
	double lMent[1];
	}TD5;
typedef TD5 **TD5Hdl;

// ============================================================================================================================================================================================================================ //
//								                                         Typedef                                                                                                                                                //
// ============================================================================================================================================================================================================================ //

// Pour les nombres de complexes
typedef std::complex<double> dcplx;

// ============================================================================================================================================================================================================================ //

namespace ROOT {
	namespace Minuit2{

enum Parameter {
				Type,
				Wavelength,
				Amplitude,
				gammaHWHM,
				dopplerHWHM,
				velocityChanging,
				parameter_0,
				parameter_1,
				parameter_2,
				parameter_3
				};


//*****************************************************************************************************************************************************************************************************************//

enum nbUseful {
			   ARRAYS_SIZE = 4,
			   NB_PARAM_LKD = 7,
			   DIM_GALATRY = 14,
			   NB_PARAM_STD = 10,
			   NB_PROFILE = 35,
			   NB_PTS_GAUSS_LEGENDRE = 512,
			   NB_PTS_GAUSS_LEGENDRE_WSL = 1024,
			   NB_PTS_GAUSS_LEGENDRE_SDG = 2048,
			   NB_PTS_GAUSS_LEGENDRE_SDG_INT = 4096			   
			  };

//*****************************************************************************************************************************************************************************************************************//

class simulation{

public:

	/**************/
	/* mutateurs  */
	/**************/

void setParam(const long &, const double &);
void setParamProfile(const std::vector<double> &,const std::vector<double> &);

	/*****************/
	/* constructeurs */
	/*****************/

simulation(const std::vector<double> &,const std::vector<double> &);
simulation(const std::vector<double> &);

/* par défaut */
simulation();

	/***************/
	/* destructeur */
	/***************/

~simulation();

/***********************/
/* Templates de calcul */
/***********************/

template <class T>
inline bool isnan(T s)
{
  // By IEEE 754 rule, 2*Inf equals Inf
  return (s!=s);
}

template <class T>
inline bool isinf(T s)
{
  // By IEEE 754 rule, 2*Inf equals Inf
  return ((2*s == s) && (s != 0));
}

// calcul du sinus cardinal 
template <class T>
inline T fcnSinc(const T x) const
{
   return (((x) != (T)0.0) ? (sin(x)/x) : ((T)1.0));
};

// calcul de x puisssance n
template <int N>
struct Power {
    template <typename NUMERIC_TYPE>
    static inline NUMERIC_TYPE of(const NUMERIC_TYPE& x) {
        return x * Power<N-1>::of(x);
    }
};
template <>
struct Power<0> {
    template <typename NUMERIC_TYPE>
    static inline NUMERIC_TYPE of(const NUMERIC_TYPE& x) {
        return 1;
    }
};

// Calcul d'une valeur inverse
template <class T>
inline T fcnInverse(T x) const
{
   return ((T)1.0/x);
};

// Calcul d'une valeur au carré
template <class T>
inline T fcnCarre(T x) const
{
   return (x*x);
};

// Calcul du théorème de Pythagore 
template <class T>
inline T fcnPythagore(const T &x,const T &y) const
{
	return (x*x + y*y);
}

// Calcul de la norme euclidienne
template <class T>
inline T fcnNorme(const T &x,const T &y) const
{
	return sqrt(x*x+y*y);
}

// Retourne le plus grand de 2 nombres
template <class T>
inline T max(const T& a,const T& b) const
{
  return (((a) >= (b)) ? (a) : (b));
};

// Retourne le plus petit de 2 nombres
template <class T>
inline T min(const T& a,const T& b) const
{
  return  (((a) <= (b)) ? (a) : (b));
};

// ============================================================================================================================================================== //

// ============================================================================================================================================================== //

	/*************************************/
	/* Fonction de simulation du spectre */
    /*************************************/

/************************************************************/
/* Fonction d'appel pour le calcul de simulation du spectre */
/************************************************************/
void profileSimulation(double *);

// Fonction de simulation pour le cas standard standard
void stdSimulationProfile(double *);

	/***********************************/
	/* SPECTRAL LINE SHAPES SIMULATION */
	/***********************************/

/******************************************/
/* BASIC MODELS FOR SPECTRAL LINE SHAPES */
/******************************************/

long doDoppler(long &, double *, std::vector<double>::const_iterator);
long doLorentz(long &, double *, std::vector<double>::const_iterator);
long doVoigt(long &, double *, std::vector<double>::const_iterator);
//long doOptimizedHumlicek(long &, double *, std::vector<double>::const_iterator);
long doHumlicek(long &, double *, std::vector<double>::const_iterator);
long doNelkinGhatak(long &, double *, std::vector<double>::const_iterator);
long doNelkinGhatakHumlicek(long &, double *, std::vector<double>::const_iterator);
long doGalatry(long &, double *, std::vector<double>::const_iterator);
long doIncompleteGamma(long &, double *, std::vector<double>::const_iterator);
long doFano(long &, double *, std::vector<double>::const_iterator);

/***********************/
/* CORRELATED PROFILES */
/***********************/

/*
long doCorrelatedNelkinGhatak(long &, double *, std::vector<double>::const_iterator);
long doPartiallyCorrelatedNelkinGhatak(long &, double *, std::vector<double>::const_iterator);
long doCorrelatedNelkinGhatakHumlicek(long &, double *, std::vector<double>::const_iterator);
long doPartiallyCorrelatedNelkinGhatakHumlicek(long &, double *, std::vector<double>::const_iterator);
long doCorrelatedGalatry(long &, double *, std::vector<double>::const_iterator);
*/

/******************************/
/* THE RAUTIAN-SOBELMAN MODEL */
/******************************/
/*
long doRautianSobelman(long &, double *, std::vector<double>::const_iterator);
long doCorrelatedRautianSobelman(long &, double *, std::vector<double>::const_iterator);
*/
/*************************************/
/* SPEED-DEPENDENT LINE-SHAPE MODELS */
/*************************************/

long doQuadraticSpeedDependentVoigt(long &, double *, std::vector<double>::const_iterator);
long doQuadraticSpeedDependentNelkinGhatak(long &, double *, std::vector<double>::const_iterator);
long doQuadraticSpeedDependentVoigtBoone(long &, double *, std::vector<double>::const_iterator);
long doQuadraticSpeedDependentNelkinGhatakBoone(long &, double *, std::vector<double>::const_iterator);
long doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone(long &, double *, std::vector<double>::const_iterator);
long doQuadraticSpeedDependentHumlicekBoone(long &, double *, std::vector<double>::const_iterator);
long doQuadraticSpeedDependentNelkinGhatakHumlicekBoone(long &, double *, std::vector<double>::const_iterator);
long doPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone(long &, double *, std::vector<double>::const_iterator);
long doApproximateQuadraticSpeedDependentGalatry(long &, double *, std::vector<double>::const_iterator);

	/************************/
	/* SPECTRAL LINE SHAPES */
	/************************/

/****************/
/* BASIC MODELS */
/****************/

double lorentzProfile(std::vector<double>::const_iterator,double &);
double dopplerProfile(std::vector<double>::const_iterator,double &);
double voigtProfile(std::vector<double>::const_iterator,double &);
//double optimizedHumlicekProfile(std::vector<double>::const_iterator,double &);
double humlicekProfile(std::vector<double>::const_iterator,double &);
double nelkinGhatakProfile(std::vector<double>::const_iterator,double &);
double nelkinGhatakHumlicekProfile(std::vector<double>::const_iterator,double &);
double galatryProfile(std::vector<double>::const_iterator,double &);
// Appendix: Computation of the Galatry Function
// [P.L. Varghese, R.K. Hanson, Collisional narrowing effects on spectral line shapes measured at high resolution.  Appl Opt 23, 2376 (1984)]
	/****************************************************************************************************************************************/
	// Pour les régions du calcul de la fonction incomplete gamma cf. P.L. Varghese and R.K. Hanson, Appl. Opt. Vol. 23, Pag. 2376, (1984). //	
	/****************************************************************************************************************************************/
	// Calcul du profil de la fonction GALATRY dans la région I
double incompleteGammaRegionI(std::vector<double>::iterator,const double &);
	// Calcul du profil de la fonction GALATRY dans la région II
double incompleteGammaRegionII(std::vector<double>::iterator,const double &);
	// Calcul du profil de la fonction GALATRY dans la région III
double incompleteGammaRegionIII(std::vector<double>::iterator,const double &);
// According to U. Fano
double fanoProfile(std::vector<double>::const_iterator,double &);

/***********************/
/* CORRELATED PROFILES */
/***********************/
/*
double correlatedNelkinGhatakProfile(std::vector<double>::const_iterator,double &);
double partiallyCorrelatedNelkinGhatakProfile(std::vector<double>::const_iterator,double &);
double correlatedNelkinGhatakHumlicekProfile(std::vector<double>::const_iterator,double &);
double partiallyCorrelatedNelkinGhatakHumlicekProfile(std::vector<double>::const_iterator,double &);
double correlatedGalatryProfile(std::vector<double>::const_iterator,double &);
*/
/******************************/
/* THE RAUTIAN-SOBELMAN MODEL */
/******************************/
/*
double rautianSobelmanProfile(std::vector<double>::const_iterator,double &);
double correlatedRautianSobelmanProfile(std::vector<double>::const_iterator,double &);
*/
/*************************************/
/* SPEED-DEPENDENT LINE-SHAPE MODELS */
/*************************************/

//-------------------------------------------------------------------------------------------
// The air-broadened, near-infrared CO2 line shape in the spectrally isolated
// regime: Evidence of simultaneous Dicke narrowing and speed dependence
// David A. Long, Katarzyna Bielska, Daniel Lisak, Daniel K. Havey, Mitchio Okumura,
// Charles E. Miller, and Joseph T. Hodges
// [THE JOURNAL OF CHEMICAL PHYSICS 135, 064308 (2011)]
//-------------------------------------------------------------------------------------------
// Parameter definitions and the quadratic approximation
// Assuming a quadratic dependence for each of the speed-
// dependent collisional functions then
// Bw(x) = 1 + aw(x*x-3/2)
// BS(x) = 1 + as(x*x-3/2)
// define the broadening speed-dependence parameter aw , and
// the shifting speed-dependence parameter as.
//-------------------------------------------------------------------------------------------
// Speed-dependent Voigt profile for water vapor in infrared remote sensing applications. 
// Boone CD, Walker KA, Bernath PF.
// J Quant Spectrosc Radiat Transf 2007;105:525-532.
//-------------------------------------------------------------------------------------------
// An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes
// N.H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann
// J Quant Spectrosc Radiat Transf 2013
//-------------------------------------------------------------------------------------------
// Efficient computation of some speed-dependent isolated line profiles
// H. Tran, N.H. Ngo, J.-M. Hartmann
// J Quant Spectrosc Radiat Transf 2013
//-------------------------------------------------------------------------------------------

double quadraticSpeedDependentVoigtProfile(std::vector<double>::const_iterator,double &);
double quadraticSpeedDependentNelkinGhatakProfile(std::vector<double>::const_iterator,double &);
double quadraticSpeedDependentVoigtBooneProfile(std::vector<double>::const_iterator,double &);
double quadraticSpeedDependentNelkinGhatakBooneProfile(std::vector<double>::const_iterator,double &);
double partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBooneProfile(std::vector<double>::const_iterator,double &);
double quadraticSpeedDependentHumlicekBooneProfile(std::vector<double>::const_iterator,double &);
double quadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(std::vector<double>::const_iterator,double &);
double partiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBooneProfile(std::vector<double>::const_iterator,double &);
double approximateQuadraticSpeedDependentGalatryProfile(std::vector<double>::const_iterator,double &);

	/********************/
	/* Calcul de l'aire */
	/********************/

/*******************************************************************/
/* Fonction d'appel pour le calcul de l'aire des différentes raies */
/*******************************************************************/

void areaProfile(double *);

/*******************************************************/
/* Fonction de calcul de l'aire des différents profils */
/*******************************************************/

/****************/
/* BASIC MODELS */
/****************/

double doAreaLorentz(std::vector<double>::const_iterator);
double doAreaDoppler(std::vector<double>::const_iterator);
double doAreaVoigt(std::vector<double>::const_iterator);
//double doAreaOptimizedHumlicek(std::vector<double>::const_iterator);
double doAreaHumlicek(std::vector<double>::const_iterator);
double doAreaNelkinGhatak(std::vector<double>::const_iterator);
double doAreaNelkinGhatakHumlicek(std::vector<double>::const_iterator);
double doAreaGalatry(std::vector<double>::const_iterator);
double doAreaIncompleteGamma(std::vector<double>::const_iterator);
double doAreaFano(std::vector<double>::const_iterator);

/***********************/
/* CORRELATED PROFILES */
/***********************/
/*
double doAreaCorrelatedNelkinGhatak(std::vector<double>::const_iterator);
double doAreaPartiallyCorrelatedNelkinGhatak(std::vector<double>::const_iterator);
double doAreaCorrelatedNelkinGhatakHumlicek(std::vector<double>::const_iterator);
double doAreaPartiallyCorrelatedNelkinGhatakHumlicek(std::vector<double>::const_iterator);
double doAreaCorrelatedGalatry(std::vector<double>::const_iterator);
*/
/******************************/
/* THE RAUTIAN-SOBELMAN MODEL */
/******************************/
/*
double doAreaRautianSobelman(std::vector<double>::const_iterator);
double doAreaCorrelatedRautianSobelman(std::vector<double>::const_iterator);
*/
/*************************************/
/* SPEED-DEPENDENT LINE-SHAPE MODELS */
/*************************************/

double doAreaQuadraticSpeedDependentVoigt(std::vector<double>::const_iterator);
double doAreaQuadraticSpeedDependentNelkinGhatak(std::vector<double>::const_iterator);
double doAreaQuadraticSpeedDependentVoigtBoone(std::vector<double>::const_iterator);
double doAreaQuadraticSpeedDependentNelkinGhatakBoone(std::vector<double>::const_iterator);
double doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakBoone(std::vector<double>::const_iterator);
double doAreaQuadraticSpeedDependentHumlicekBoone(std::vector<double>::const_iterator);
double doAreaQuadraticSpeedDependentNelkinGhatakHumlicekBoone(std::vector<double>::const_iterator);
double doAreaPartiallyCorrelatedQuadraticSpeedDependentNelkinGhatakHumlicekBoone(std::vector<double>::const_iterator);
double doAreaApproximateQuadraticSpeedDependentGalatry(std::vector<double>::const_iterator);


// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

	/********************************/
	/* Calcul pour la ligne de base */
	/********************************/

/*******************************************************/
/* Fonction d'appel pour le calcul de la ligne de base */
/*******************************************************/
void FCN_Baseline(double *,const std::vector<double> &);
inline void BaselineConstant(double *, const std::vector<double> &);
void BaselinePolynomial1(double *, const std::vector<double> &, std::vector<double>::const_iterator);
void BaselinePolynomial2(double *, const std::vector<double> &, std::vector<double>::const_iterator);
void BaselinePolynomial3(double *, const std::vector<double> &, std::vector<double>::const_iterator);
void BaselinePolynomial4(double *, const std::vector<double> &, std::vector<double>::const_iterator);
void BaselinePolynomial5(double *, const std::vector<double> &, std::vector<double>::const_iterator);
void BaselinePolynomialn(double *, const std::vector<double> &, std::vector<double>::const_iterator);

	// Fonction pour calculer les abscisses pour la ligne de base compris dans [-1;1]
inline void initAbscisse(std::vector<double> &);

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

																				// ********************************************** //
																				/*         Fonctions auxiliaires                  */
																				// ********************************************** //

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

  // ======================================================================================== //
  //                     Fonction erreur complexe (CERN)
  // ======================================================================================== //
  /*                                                                                          */
  /* Return CERNlib complex error function                                                    */
  /*                                                                                          */
  // ======================================================================================== //

dcplx complexErrFunc(const dcplx &);

  // ======================================================================================== //
  /*                                                                                          */
  /* This Routine was Taken from the Paper by J. Humlicek, which                              */
  /* is Available in Page 309 of Volume 21 of the 1979 Issue of                               */
  /* the Journal of Quantitative Spectroscopy and Radiative Transfer                          */
  /*                                                                                          */
  // ======================================================================================== //

dcplx CPF(const dcplx &);
dcplx CPF3(const dcplx &);

	/**********************************************************************************************/
	//           J. Humlicek, J Quant. Spectrosc. Radiat. Transfer, vol. 21,pag. 309 (1979)       //
	/**********************************************************************************************/
	//  Calcul, pour la position z=x+iy, de la fonction de VOIGT suivant l'algorithme d'HUMLICEK  //
	/**********************************************************************************************/

//dcplx humlicekValue(const dcplx &);

	/**********************************************************************************************/
	//            Approximations utilisées pour le calcul de la fonction HUMLICEKVALUE            //
	/**********************************************************************************************/

//dcplx approximation1(const dcplx &);
//dcplx approximation2(const dcplx &,const dcplx &);
//dcplx approximation3(const dcplx &);
//dcplx approximation4(const dcplx &,const dcplx &);

// ======================================================================================== //
// Calcul des coefficients de l'expansion asymptotique pour la fonction GALATRY 
// ======================================================================================== //

void coeffAsymtoticExpansion(std::vector<double>::iterator);

// ======================================================================================== //
// Initialisation des valeurs pour le calcul de l'intégrale dans le cas 
// "quadratic speed-dependent"
// ======================================================================================== //

void initQSDVP(const long &,std::vector<double>::iterator);
void initQSDNGP(const long &,std::vector<double>::iterator);
void initQSDGP(const long &,std::vector<double>::iterator);

// ============================================================================================================================================================================================================================ //
// ============================================================================================================================================================================================================================ //

private:
// données membres
	std::vector<double> fPar;
	std::vector<double> fX;
};
	}  // namespace Minuit2
} // namespace ROOT
#endif // MN_SIMULATION_H