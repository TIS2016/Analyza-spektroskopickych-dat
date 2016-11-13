/* Call Library source file */

#include "extcode.h"
#include "Windows.h"
#include <string>
#include <ppl.h>
/* lv_prolog.h and lv_epilog.h set up the correct alignment for LabVIEW data. */
#include "lv_prolog.h"

using namespace concurrency;
using namespace std;

/* Typedefs */

typedef struct {
	int32_t dimSizes[2];
	LStrHandle String[1];
} TD2;
typedef TD2 **TD2Hdl;

typedef struct {
	int32_t dimSizes[2];
	double Numeric[1];
} TD3;
typedef TD3 **TD3Hdl;

typedef struct {
	int32_t dimSizes[3];
	int32_t Numeric[1];
} TD4;
typedef TD4 **TD4Hdl;

typedef struct {
	int32_t dimSize;
	LStrHandle String[1];
} TD5;
typedef TD5 **TD5Hdl;

typedef struct {
	int32_t dimSize;
	int32_t Numeric[1];
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
	int32_t dimSize;
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
	int32_t code;
	LStrHandle source;
} TD10;

#include "lv_epilog.h"

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TD9 *DataOUT, TD10 *Error);

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TD9 *DataOUT, TD10 *Error)
{

	/* Insert code here */

	return 7;
}

