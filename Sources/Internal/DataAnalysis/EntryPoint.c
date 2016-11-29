#include "stdafx.h"
#include "EntryPoint.h"

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TDFast *DataOUT_F, TD10 *Error) {
	int parameterCount = 0;
	{
		DebugLogger logger("F:\\School\\3.rocnik\\Winter\\TIS\\Log01.txt");
		logger.Log("logger OK!");

		shared_ptr<TD1> spInputParams(DataPARin);
		logger.Log("smartPointer OK!");

		parameterCount = 9;//*((*(DataPARin->Data_length))->Numeric);
		logger.Log("paramCount OK?");

		string tmp;
		LStrToStr(*spInputParams->Name, tmp);
		logger.Log("string conversion OK!");

	}
	UNREFERENCED_PARAMETER( DataPARin );
	UNREFERENCED_PARAMETER(DataIN);
	UNREFERENCED_PARAMETER(DataOUT_F);
	UNREFERENCED_PARAMETER(Error);
	return parameterCount + 5;
}

extern "C" _declspec(dllexport) int32_t fdata_complete(TD1 *DataPARin, TD7 *DataIN, TDComplete *DataOUT_C, TD10 *Error) {
	UNREFERENCED_PARAMETER(DataPARin);
	UNREFERENCED_PARAMETER(DataIN);
	UNREFERENCED_PARAMETER(DataOUT_C);
	UNREFERENCED_PARAMETER(Error);
	return 88;
}