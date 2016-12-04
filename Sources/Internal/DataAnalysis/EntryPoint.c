#include "stdafx.h"
#include "EntryPoint.h"
#include <vector>

void DumpToLog(__in vector<string> &toDump, __in shared_ptr<DebugLogger> spLogger) { 
	spLogger->BeginSection();
	spLogger->Log("Dumping...");
	for (auto it = toDump.begin(); it != toDump.end(); it++) {
		spLogger->Log(*it);
	}
	spLogger->EndSection();
}

void GetInputParamStrings( __in shared_ptr<TD1> spInParams, __out vector<string> &outParams) {
	string tmp;
	TD2Hdl firstParamStrDeref = spInParams->ParamStrings;
	size_t count = (*firstParamStrDeref)->dimSizes[0];

	TD2 *pParamIn = *firstParamStrDeref;
	for (size_t i = 0; i < count; i++) {
		LStrPtr inStr = **(pParamIn->String);
		LStrToStr(inStr, tmp);

		// TODO: trim
		outParams.push_back(tmp);

		tmp.clear();
		pParamIn++;
	}
}

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TDFast *DataOUT_F, TD10 *Error) {
	int parameterCount = 0;
	{
		shared_ptr<DebugLogger> spLogger( new DebugLogger(LOG_PATH_FILIP) );

		spLogger->Log("logger OK!");

		shared_ptr<TD1> spInputParams(DataPARin);
		spLogger->Log("smartPointer OK!");

		string dumpStr;
		LStrToStr(*(spInputParams->Name), dumpStr);
		spLogger->Log(dumpStr);

		if ( spInputParams->Data_length != nullptr) {
			spLogger->LogFormatted("Data_length: %d", 9 );
		} else {
			spLogger->Log("Data_length is nullptr");
		}

		vector<string> paramStrings;
		paramStrings.reserve(256);
		GetInputParamStrings(spInputParams, paramStrings);
		DumpToLog(paramStrings, spLogger);

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