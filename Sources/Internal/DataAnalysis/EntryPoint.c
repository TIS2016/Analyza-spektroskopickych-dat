/* Call Library source file */

#include "stdafx.h"
/* lv_prolog.h and lv_epilog.h set up the correct alignment for LabVIEW data. */
#include "../../External/cintools/lv_prolog.h"

/* Typedefs */

using namespace std;
using namespace concurrency;

typedef struct {
	INT dimSizes[2];
	LStrHandle String[1];
} TD2;
typedef TD2 **TD2Hdl;

typedef struct {
	INT dimSizes[2];
	double Numeric[1];
} TD3;
typedef TD3 **TD3Hdl;

typedef struct {
	INT dimSizes[3];
	INT Numeric[1];
} TD4;
typedef TD4 **TD4Hdl;

typedef struct {
	INT dimSize;
	LStrHandle String[1];
} TD5;
typedef TD5 **TD5Hdl;

typedef struct {
	INT dimSize;
	INT Numeric[1];
} TD6;
typedef TD6 **TD6Hdl;

typedef struct {
	LStrHandle Name;
	TD2Hdl ParamStrings;
	TD3Hdl ParamNumbers;
	TD2Hdl FuncNames;
	TD4Hdl FuncParAdresses;
	TD5Hdl DataNames;
	TD6Hdl Data_length;
} TD1;

typedef struct {
	INT dimSize;
	double Numeric[1];
} TD8;
typedef TD8 **TD8Hdl;

typedef struct {
	LStrHandle Name;
	TD8Hdl XIn;
	TD8Hdl YIn;
	TD8Hdl WIn;
} TD7;

typedef struct {
	LStrHandle Name;
	TD8Hdl XOut;
	TD8Hdl YOut;
	TD8Hdl WOut;
	TD8Hdl F;
	TD3Hdl _2DData;
} TD9;

typedef struct {
	LVBoolean status;
	INT code;
	LStrHandle source;
} TD10;

#include "../../External/cintools/lv_epilog.h"

// extern "C" _declspec(dllexport) INT fdata_fast(TD1 *DataPARin, TD7 *DataIN, TD9 *DataOUT, TD10 *Error);

extern "C" _declspec(dllexport) INT fdata_fast(TD1 *DataPARin, TD7 *DataIN, TD9 *DataOUT, TD10 *Error)
{
	UNREFERENCED_PARAMETER( DataPARin );
	UNREFERENCED_PARAMETER(DataIN);
	UNREFERENCED_PARAMETER(DataOUT);
	UNREFERENCED_PARAMETER(Error);
	return 7;
}

