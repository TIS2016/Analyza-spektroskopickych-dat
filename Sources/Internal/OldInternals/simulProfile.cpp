#include "simulation.h"
#include "extcode.h"
#include <ppl.h>
#include <vector>
#include "concurrent_vector.h"

using namespace std;
using namespace Concurrency;

using namespace ROOT::Minuit2;




//#include "C:\Users\Serge\Documents\Visual Studio 2010\LIBRARY\CINTOOLS_LV_7_1\extcode.h"


/* Typedefs */

//Peter
typedef struct {
	LVBoolean AmplitudeDriven;
	LVBoolean CalculateParameterWeight;
} pTD1;

typedef struct {
	long dimSize;
	double elt[1];
} pTD2;
typedef pTD2 **pTD2Hdl;

typedef struct {
	long dimSizes[2];
	double elt[1];
} pTD3;
typedef pTD3 **pTD3Hdl;


// Serge
typedef struct {
	double Wavelength;
	double Amplitude;
	double Lorentz;
	double Gauss;
	double Gamma;
	double param0;
	double param1;
	double param2;
	double param3;
	long FitType;
	LVBoolean WLFit;
	LVBoolean AmpFit;
	LVBoolean LorFit;
	LVBoolean GaussFit;
	LVBoolean GamFit;
	LVBoolean param0Fit;
	LVBoolean param1Fit;
	LVBoolean param2Fit;
	LVBoolean param3Fit;
	} TD2;

typedef struct {
	long dimSize;
	TD2 Cluster2[1];
	} TD1;
typedef TD1 **TD1Hdl;

// Baseline parameters
typedef struct {
	double a;
	LVBoolean fit;
	} TD4;

typedef struct {
	long dimSize;
	TD4 Cluster[1];
	} TD3;
typedef TD3 **TD3Hdl;




/**********************************************************************************************************************/
/*											LINE SIMULATION  Peter  		            									  */
/**********************************************************************************************************************/

/* For Parallel case...*/

_declspec(dllexport) int32 Peaks
(
	pTD1 *Options,
	pTD2Hdl *X,
	pTD2Hdl *Y,
	pTD3Hdl *PeaksData,
	pTD3Hdl *PeaksDataWeight
	)

{
	/**************************************************/
	// 1 = "No data"
	// 2 = "Bad Y initialisation"
	// 3 = "No peaks"
	// 4 = "LinesM error"
	// Checking

	size_t PointsNumber, PeaksNumber, ParPerPeak;

	double *PointerToX, *PointerToY;

	if ((*X) == NULL) return 1;
	if ((*Y) == NULL) return 2;

	PointsNumber = (**X)->dimSize;
	PeaksNumber = (**PeaksData)->dimSizes[0];     //this is Serges'nblines
	ParPerPeak = (**PeaksData)->dimSizes[1];

	PointerToX = (**X)->elt;
	PointerToY = (**Y)->elt;

	if (PeaksNumber == 0) return 3;

	//put the data in vector container
	//std::vector<double> XscaleVec(PointerToX, PointerToX + PointsNumber);
	//std::vector<double> YscaleVec(PointerToY, PointerToY + PointsNumber);
	//vector<double> PeaksParam(PeaksData, PeaksNumber*ParPerPeak );
	//std::vector<double>::iterator yptr;
	//std::vector<double>::iterator xptr;

	//yptr = YscaleVec.begin();
	//xptr = XscaleVec.begin();

	for (size_t i = 0; i < PointsNumber; i++)
	{
		for (size_t j = 0; j < 20; j++)
		{
			PointerToY[i] += 0.001 / (0.001 + (PointerToX[i] - (6000.5 + j*0.1))*(PointerToX[i] - (6000.5 + j*0.1)));
		};
	};


//	yptr = YscaleVec.begin(); //point back to first element

//	for (size_t i = 0; i < PointsNumber; i++)
//	{
//		PointerToY[i] = *yptr;
//		++yptr;
//	}




	return 1000;
}


_declspec(dllexport) int32 PeaksPar 
(
	pTD1 *Options,	
	pTD2Hdl *X,
	pTD2Hdl *Y,
	pTD3Hdl *PeaksData,
	pTD3Hdl *PeaksDataWeight
	)

{
	
	/**************************************************/
	// 1 = "No data"
	// 2 = "Bad Y initialisation"
	// 3 = "No peaks"
	// 4 = "LinesM error"
	// Checking

	size_t PointsNumber, PeaksNumber, ParPerPeak;

	double *PointerToX, *PointerToY;

	if ((*X) == NULL) return 1;
	if ((*Y) == NULL) return 2;

	PointsNumber = (**X)->dimSize;
	PeaksNumber = (**PeaksData)->dimSizes[0]; //this is Serges'nblines
	ParPerPeak = (**PeaksData)->dimSizes[1];

	PointerToX = (**X)->elt;
	PointerToY = (**Y)->elt;

	if (PeaksNumber == 0) return 3;

	//double dummy ()

	//put the data in vector container
	//concurrent_vector<double> XscaleVec(PointerToX, PointerToX+PointsNumber);
	//concurrent_vector<double> YscaleVec(PointerToY, PointerToY+PointsNumber);
	//vector<double> PeaksParam(PeaksData, PeaksNumber*ParPerPeak );
	//concurrent_vector<double>::iterator yptr;
	//concurrent_vector<double>::iterator xptr;


	parallel_for(size_t(0), PointsNumber, [&](size_t i) 
	{
		for (size_t j = 0; j < 20; j++) 
		{
		PointerToY[i] += 0.001 / (0.001 + (PointerToX[i] - (6000.5+j*0.1))*(PointerToX[i] - (6000.5 + j*0.1)));
		};
		//YscaleVec[i] = 0.5 / (0.5 + (XscaleVec[i] - 6002.5)*(XscaleVec[i] - 6002.5));
	}
	);
	
	//now write the result into LB array:)

	//parallel_for_each(YscaleVec.begin(), YscaleVec.end(), [&](double y)
	//{
	//});
	
//	yptr = YscaleVec.begin();
//	xptr = XscaleVec.begin();
	
//	for (size_t i=0; i < PointsNumber; i++) {
//		*yptr = 0.5/(0.5+(*xptr-6002.5)*(*xptr - 6002.5));
//		++xptr;
//		++yptr;
//	}



//	yptr = YscaleVec.begin(); //point back to first element

//	for (size_t i=0; i < PointsNumber; i++)
//	{
//		PointerToY[i] = *yptr;
//		++yptr;
//	}
	
	return 0;
}


/**********************************************************************************************************************/
/*											LINE SIMULATION Serge	            									  */
/**********************************************************************************************************************/


__declspec(dllexport) int32 _lineSimulation (long nbData, double xData[], double simulData[], double simulBaseline[],
										   double simulDataWithoutBaseline[], TD1Hdl *lineList, TD3Hdl *Baseline,
										   unsigned short typeA, double area[])
{

	long nbLines(0);
	long nbCoeffsPolynome(0);

	/*------------------------------------------------------------------------------*/
	// Si le tableau de cluster est vide (càd (**X)->dimSize==0), le ptr est égal à 0:
	// LABVIEW ne crée pas une structure de dimension 0 pour un
	// passage par pointeur vers handles contrairement au passage
	// handles par valeur.
	// Donc l'opération de lecture (**Baseline)->dimSize est impossible
	// car le cluster n'a pas été crée
	/*------------------------------------------------------------------------------*/

	if (*lineList != NULL && *Baseline != NULL)
	{
		nbLines = (**lineList)->dimSize;
		nbCoeffsPolynome=(**Baseline)->dimSize;
	}
	else if (*lineList != NULL && *Baseline == NULL)
		nbLines = (**lineList)->dimSize;
	else if (*lineList == NULL && *Baseline != NULL)
		nbCoeffsPolynome=(**Baseline)->dimSize;
	else
		return -1;

	long nbParametres = ARRAYS_SIZE + NB_PARAM_STD*nbLines + nbCoeffsPolynome;
	double gHWHM(0.0);
	vector<double> abscisses(xData,xData+nbData);
	vector<double> par(nbParametres,0.);
	vector<double>::iterator ptrParam = par.begin()+ARRAYS_SIZE;


	par[0] = (double)nbLines;				// #raies
	par[1] = (double)nbCoeffsPolynome;		// #coefficients de la ligne de base
	par[2] = (double)typeA;					// Amplitude or area
	par[3] = (double)0.;					// Standard method
	
	for (long i=0;i<nbLines;++i)
	{
		*ptrParam = ((double)(**lineList)->Cluster2[i].FitType);
		++ptrParam;
		*ptrParam = (**lineList)->Cluster2[i].Wavelength;
		++ptrParam;
		*ptrParam = (**lineList)->Cluster2[i].Amplitude;
		++ptrParam;
		*ptrParam = (**lineList)->Cluster2[i].Lorentz;
		++ptrParam;
		gHWHM = ((**lineList)->Cluster2[i].Gauss)*INV_SQRT_LN2;
		*ptrParam = gHWHM;
		++ptrParam;
		*ptrParam = (**lineList)->Cluster2[i].Gamma;
		++ptrParam;
		*ptrParam = ((double)(**lineList)->Cluster2[i].param0);
		++ptrParam;
		*ptrParam = ((double)(**lineList)->Cluster2[i].param1);
		++ptrParam;
		*ptrParam = ((double)(**lineList)->Cluster2[i].param2);
		++ptrParam;
		*ptrParam = ((double)(**lineList)->Cluster2[i].param3);
		++ptrParam;
	}

	if (nbCoeffsPolynome)
	{
		for (int i=0;i<nbCoeffsPolynome;++i)
		{
			*ptrParam = (**Baseline)->Cluster[i].a;
			++ptrParam;
		}
	}

	//*******************//
	// Calcul simulation //
	//*******************//

	simulation fcn(abscisses,par);

	// Simulation des data
	fcn.profileSimulation(simulData);
    // Simulation de la ligne de base
	if (nbCoeffsPolynome || (nbCoeffsPolynome == 1 && (**Baseline)->Cluster[0].a != (double)0.0))
		fcn.FCN_Baseline(simulBaseline,par);
	// Simulation du spectre sans la ligne de base
	// mise à 0 du nombre de coefficients du polynôme
	// de la ligne de base
	fcn.setParam(1,0);
	// Simulation du spectre
	fcn.profileSimulation(simulDataWithoutBaseline);

	if (!typeA)
		fcn.areaProfile(area);

	return 0;

}
