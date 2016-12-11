#include "stdafx.h"
#include "EntryPoint.h"
#include <vector>

struct Sample {
	double X;
	double Y;
	double Deviation;
};

struct Function {
	string name;
	vector<string> subFunctions;
};

template <typename T> __forceinline void FillArr(__in size_t count, __in_ecount(count) T *pSrc, __out_ecount(count) T *pDst) {
	memcpy(pDst, pSrc, count * sizeof(T));
}

template <typename LWStruct> __forceinline LWStruct* LwArrayGet(__in LWStruct *pSrc, __in size_t rowCount, __in size_t row, __in size_t col) {
	return ( pSrc + col * rowCount + row );
}

template <typename LWStruct> __forceinline LWStruct* LwArrayGet(__in LWStruct *pSrc, __in size_t index ) {
	return ( pSrc + index );
}

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

__forceinline void GetInputParameterStrings(__in size_t stringCount, __in_ecount(stringCount) TD2 *pInputStrings, __out_ecount(stringCount) string *pOut) {
	for (size_t i = 0; i < stringCount; i++) {
		LStrPtr inStr = **(pInputStrings->String);
		LStrToStr(inStr, *pOut);

		pOut++;
		pInputStrings++;
	}
}

__forceinline void GetInputParameterValues(__in size_t paramCount, __in_ecount(paramCount) TD3 *pSrc, __out_ecount(paramCount) double *pDst) {
	memcpy( pDst, pSrc->Numeric, paramCount * sizeof(double) );
}

void ConvertInputData(__in size_t dataCount, __in_ecount(dataCount) double *pXin,
	__in_ecount(dataCount) double *pYin, __in_ecount(dataCount) double *pDeviationIn,
	__out_ecount(dataCount) Sample *pOut) {
	for (size_t i = 0; i < dataCount; i++) {
		pOut->X = *pXin;
		pOut->Y = *pYin;
		pOut->Deviation = *pDeviationIn;

		pXin++;
		pYin++;
		pDeviationIn++;
		pOut++;
	}
}

void GetSubFunctions(string &pSrc, vector<string> &pDst) {
	auto begin = pSrc.begin();
	size_t lastPos = 0;

	while ( true ) {
		lastPos = pSrc.find('@', lastPos);
		if (lastPos == string::npos) {
			auto end = pSrc.end();
			pDst.push_back(string(begin, end));
			break;
		} else {
			auto end = pSrc.begin() + lastPos;
			pDst.push_back(string(begin, end));
			begin = end + 1;
		}
		lastPos++;
	}
}

void GetFunctionsToApply(__in size_t count, __in_ecount(count) TD2 *pSrcStruct, __out_ecount(count) Function *pDst) {
	LStrHandle *pSrc = pSrcStruct->String;

	for (size_t i = 0; i < count; i++) {

		LStrToStr(**(pSrc), (pDst->name));
		
		string subFuncts;
		LStrToStr(**(pSrc + 1), subFuncts);
		LTrim(subFuncts, " \t\n\v\f\r");
		RTrim(subFuncts, " \t\n\v\f\r");

		GetSubFunctions(subFuncts, pDst->subFunctions);
		
		LTrim(pDst->name, " \t\n\v\f\r");
		RTrim(pDst->name, " \t\n\v\f\r");

		pSrc += 2;
		pDst++;
	}
}

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TDFast *DataOUT_F, TD10 *Error) {
	TD1 *pInputParametes = DataPARin;
	UNREFERENCED_PARAMETER(pInputParametes);
	TD7 *pInputData = DataIN;
	UNREFERENCED_PARAMETER(pInputData);
	TDFast *pOutputData = DataOUT_F;
	UNREFERENCED_PARAMETER(pOutputData);

	size_t parameterCount = static_cast<size_t>( (*(pInputParametes->ParamStrings))->dimSizes[0] );

	Buffer<string> parameterStrings;
	parameterStrings.Allocate(parameterCount);
	GetInputParameterStrings( parameterCount, *(pInputParametes->ParamStrings), parameterStrings.Ptr() );

	Buffer<double> parameterValues;
	parameterValues.Allocate(parameterCount);
	GetInputParameterValues( parameterCount, *(pInputParametes->ParamNumbers), parameterValues.Ptr() );

	size_t inputSamplesCount = static_cast<size_t>( (*(pInputData->X_in))->dimSize );
	
	Buffer<Sample> inputSamples;
	inputSamples.Allocate(inputSamplesCount);
	ConvertInputData(
		inputSamplesCount,
		(*(pInputData->X_in))->Numeric,
		(*(pInputData->Y_in))->Numeric,
		(*(pInputData->W_in))->Numeric,
		inputSamples.Ptr()
	);

	size_t functionCount = static_cast<size_t>( (*(pInputParametes->FuncNames))->dimSizes[0] );
	Buffer<Function> functsToApply;
	functsToApply.Allocate(functionCount);
	GetFunctionsToApply(functionCount, *pInputParametes->FuncNames, functsToApply.Ptr());
	
	TD4 *functParamAdresses = *(pInputParametes->FuncParAdresses);
	TD3 *paramValues = *(pInputParametes->ParamNumbers);

	size_t functionParamAdressSpace = functParamAdresses->dimSizes[1] * functParamAdresses->dimSizes[2];
	size_t subGroupAddressSpace = functParamAdresses->dimSizes[2];

	INT *pAddresses = functParamAdresses->Numeric;

	struct AddrValPair {
		INT addr;
		DOUBLE val;
	};
	vector<AddrValPair> addresses;
	
	for (size_t functI = 0; functI < functionCount; functI++) {
		Function funct = functsToApply[functI];
		INT *pAddr = pAddresses;
		for (size_t functGroupI = 0; functGroupI < funct.subFunctions.size(); functGroupI++) {
			
			size_t paramAddrCount = *pAddr;
			for (size_t paramAddrI = 0; paramAddrI < paramAddrCount; paramAddrI++) {
				INT addr = *(pAddr + 1 + paramAddrI);
				DOUBLE val = *(paramValues->Numeric + 3 * addr );
				addresses.push_back( AddrValPair{ addr, val } );
			}

			pAddr += subGroupAddressSpace;
		}

		pAddresses += functionParamAdressSpace;
	}
	

	UNREFERENCED_PARAMETER(Error);
	return 343;
}

extern "C" _declspec(dllexport) int32_t fdata_complete(TD1 *DataPARin, TD7 *DataIN, TDComplete *DataOUT_C, TD10 *Error) {
	UNREFERENCED_PARAMETER(DataPARin);
	UNREFERENCED_PARAMETER(DataIN);
	UNREFERENCED_PARAMETER(DataOUT_C);
	UNREFERENCED_PARAMETER(Error);
	return 88;
}