/* Call Library source file */
/* developped by Peter Cermak for the spectra synthetizator
 * this is the version with simple arrays input/output

 */
#include "Windows.h"
#include "extcode.h"
#include <string>
#include <ppl.h>
//#include <cmath>
//#include <complex>

using namespace concurrency;
using namespace std;

//#include <stdlib.h>
//#include <stdio.h>

/* lv_prolog.h and lv_epilog.h set up the correct alignment for LabVIEW data. */
#include "lv_prolog.h"

/* Typedefs */

typedef struct {
	LVBoolean AmplitudeDriven;
	LVBoolean CalculateParameterWeight;
	} TD1;

typedef struct {
	long dimSize;
	double elt[1];
	} TD2;
typedef TD2 **TD2Hdl;

typedef struct {
	long dimSizes[2];
	double elt[1];
	} TD3;
typedef TD3 **TD3Hdl;

typedef struct {
	long dimSizes[2];
	int32 elt[1];
} TD4;
typedef TD4 **TD4Hdl;

#include "lv_epilog.h"
#include "linesAdv.h"
#include "Cerf.h"

/* This is the 2016 Function...*/

extern "C" _declspec(dllexport) int32_t Peaks_2016
(
	int32_t Options,
	TD2Hdl *X,
	TD2Hdl *Y,
	TD3Hdl *PeaksData,
	TD3Hdl *PeaksAdr
	)

{



	return 0;
}


/* This is the Parallel Function...*/

extern "C" _declspec(dllexport) int32_t Wcalc
(
	int32_t Options,
	TD2Hdl *X,
	TD2Hdl *Y
	)

{
	size_t i, PointsNumber;
	double vr(0.0);
	double vi(0.0);
		
	double *PointerToX, *PointerToY;

	PointsNumber = (**X)->dimSize;
	PointerToX = (**X)->elt;
	PointerToY = (**Y)->elt;


	switch (Options)
	{
			//Marco's CERN representation***********************
		case 0: {
			for (i = 0; i < PointsNumber; i++)
					{
				cwerf(PointerToX[i], PointerToY[i], &vr, &vi);
				PointerToX[i] = vr;
				PointerToY[i] = vi; // vi;
					}
				}
			break; 

			//Serge's CERN representation***********************
		case 1: {
			for (i = 0; i < PointsNumber; i++)
					{
				dcplx z(PointerToX[i], PointerToY[i]);
				dcplx w(complexErrFunc(z));
				PointerToX[i] = w.real();
				PointerToY[i] = w.imag();
					}
				}
				break;
			
				//Serge's HUM representation***********************
		case 2: {
			for (i = 0; i < PointsNumber; i++)
					{
				dcplx z(PointerToX[i], PointerToY[i]);
				dcplx w(CPF(z));
				PointerToX[i] = w.real();
				PointerToY[i] = w.imag();
					}
				}
			break;

			//CERN's HUM representation***********************
		case 3: {
			for (i = 0; i < PointsNumber; i++)
			{
				dcplx z(PointerToX[i], PointerToY[i]);
				dcplx w(Cerf::faddeeva(z));
				PointerToX[i] = w.real();
				PointerToY[i] = w.imag();
			}
		}
				break;
			
		default:
			break;
	};

	return 0;
}

/* This is the Parallel Function...*/

extern "C" _declspec(dllexport) int32_t PeaksPar
(
	int32_t Options,
	TD2Hdl *X,
	TD2Hdl *Y,
	TD3Hdl *PeaksData
	)


{
	/**************************************************/
	/* Local variables                                */

	int shape;
	bool am;
	size_t PointsNumber, PeaksNumber, ParNumber;
//	double buffer1, buffer2;
	double *PointerToX, *PointerToY, *PointerToPeak;

	/*Shortcuts for some variables*********************/
	
	PointsNumber = (**X)->dimSize;
	PeaksNumber = (**PeaksData)->dimSizes[0];
	ParNumber = (**PeaksData)->dimSizes[1];
	
	if (Options == 0) am = true;
	else am = false;
	/**************************************************/
	// 1 = "No data"
	// 2 = "Bad Y initialisation"
	// 3 = "No peaks"
	// 4 = "LinesM error"
	// Checking

	/*Check for errors*********************************/
	if ((*X) == NULL) return 1;
	if ((*Y) == NULL) return 2;
    if (PeaksNumber == 0) return 3;
	
	/* Setting pointers********************************/
	PointerToX = (**X)->elt;
	PointerToY = (**Y)->elt;
	PointerToPeak = (**PeaksData)->elt;

	/* Point by point calculation */
	parallel_for(size_t(0), PointsNumber, [=](size_t i)
	{
		for (size_t p=0; p < PeaksNumber; p++) 
		{ 
			PointerToY[i] += PeakAdv(PointerToX[i], PointerToPeak+p*ParNumber, am);
		};    
	});

	/*end of simulation*/
	return 0;
}





/* This is the Serial Function...*/

extern "C" _declspec(dllexport) int32_t Peaks
(
	int32_t Options,
	TD2Hdl *X, 
	TD2Hdl *Y, 
	TD3Hdl *PeaksData
	
)

{
	/**************************************************/
	// Local variables 

	int shape;
	bool am;
	size_t PointsNumber, PeaksNumber, ParNumber;
//	double buffer1, buffer2;
	double *PointerToX, *PointerToY, *PointerToPeak;
    
	
	/*Shortcuts for some variables*/

	PointsNumber = (**X)->dimSize;
	PeaksNumber = (**PeaksData)->dimSizes[0];
	ParNumber = (**PeaksData)->dimSizes[1];

	if (Options == 0) am = true;
	else am = false;
	
	/*Check for errors**********************************************/
	// 1 = "No data"
	// 2 = "Bad Y initialisation"
	// 3 = "No peaks"
	// 4 = "LinesM error"
	// Checking

	if ((*X) == NULL) return 1;
	if ((*Y) == NULL) return 2;
	if (PeaksNumber == 0) return 3;
		

	/* Setting pointers */
	PointerToX = (**X)->elt;
	PointerToY = (**Y)->elt;
	PointerToPeak = (**PeaksData)->elt;
	
	/* Point by point calculation */
	for (size_t i = 0; i < PointsNumber; i++)
	{
		for (size_t p = 0; p < PeaksNumber; p++)
		{
			PointerToY[i] += PeakAdv(PointerToX[i], PointerToPeak + p*ParNumber, am);
		};
	}
	
	/*/Add each peak
	for (i=0; i<PeaksNumber; i++)
	{
		// in PeaksData each peak corresponds to one line containing type(0), center(1), ampl.(2), wG(3), wL(4), wC(5), surface(6) (7 parameters)
		PointerToPeak = ((**PeaksData)->elt)+i*7;
		
		// to widths used in LinesM:
		shape = int(*PointerToPeak);
		buffer1 = *(PointerToPeak+3)/(2*sqrt(log(2)));  //from FWHM to s...: s = FWHM/(2 sqrt(ln(2)))
		buffer2 = *(PointerToPeak+4)/2;               //from FWHM to g...: g = 1/2 FWHM

		switch(shape)
			{
				//lorentz***********************
			case 0: 
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / lornorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * lornorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( lorentz( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

				//gauss*************************
			case 1:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / gaunorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * gaunorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( gauss( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;
			
				//voigt*************************
			case 2:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / voinorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * voinorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( voigt( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

				//humlicek*************************
			case 3:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / voinorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * voinorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( humlicek( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

				//rautian*************************
			case 4:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / raunorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * raunorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( rautian( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

				//rautianh*************************
			case 5:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / raunorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * raunorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( rautianh( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

				//galatry*************************
			case 6:
				//calc. peak amplitude if not ampl. driven: ampl= surf./lornorm
				if (Options->AmplitudeDriven != 1)  *(PointerToPeak+2) = *(PointerToPeak+6) / galnorm(buffer2, buffer1, *(PointerToPeak+5));
					else *(PointerToPeak+6) = *(PointerToPeak+2) * galnorm(buffer2, buffer1, *(PointerToPeak+5));
				//peak
				if ( galatry( PointsNumber, PointerToX, PointerToY, *(PointerToPeak+1), *(PointerToPeak+2), buffer2, buffer1, *(PointerToPeak+5)) != 0) return 4;
		   
				break;

			}
	
	} */
	
	return *(PointerToPeak +0);
}
