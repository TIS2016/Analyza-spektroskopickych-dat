#pragma once

/* Call Library source file */

/* lv_prolog.h and lv_epilog.h set up the correct alignment for LabVIEW data. */
#include "../../External/cintools/lv_prolog.h"

/* Typedefs */

using namespace std;
using namespace concurrency;

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
	TD8Hdl X_in;
	TD8Hdl Y_in;
	TD8Hdl W_in;
} TD7;

typedef struct {
	LStrHandle Name;
	TD8Hdl X_out;
	TD8Hdl Y_out;
	TD8Hdl W_out;
	TD8Hdl F;
} TDFast;

typedef struct {
	LStrHandle Name;
	TD8Hdl X_out;
	TD8Hdl Y_out;
	TD8Hdl W_out;
	TD8Hdl F;
	TD5Hdl _2DData_In_Names;
	TD3Hdl _2DData_In;
	TD5Hdl _2DData_Out_Names;
	TD3Hdl _2DData_Out;
} TDComplete;

typedef struct {
	LVBoolean status;
	int32_t code;
	LStrHandle source;
} TD10;

#include "../../External/cintools/lv_epilog.h"

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TDFast *DataOUT_F, TD10 *Error);

extern "C" _declspec(dllexport) int32_t fdata_complete(TD1 *DataPARin, TD7 *DataIN, TDComplete *DataOUT_C, TD10 *Error);