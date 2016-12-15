#pragma once

#include "EntryPoint.h"

namespace DataAnalysis { namespace InputOutput {

	inline void GetTransformationName( __in LStrHandle *pName, __out FunctionHeader *pHeader );

	inline void SplitSubfunctions( __in string &subFunctString, __out FunctionHeader *pHeader );

	inline void GetTransformationSubfunctions( __in LStrHandle *pString, __out FunctionHeader *pHeader );

	void GetTransformationHeaders( __in size_t count, __in_ecount( count ) TD2 *pSrcStruct, __out_ecount( count ) FunctionHeader *pDst );

	void GetFunctionParameters( __in TD4 *pParameterIndexes, __in TD3 *pParameters, __inout Buffer<FunctionHeader> &headers );

	void ConvertInputData( __in TD7 *pInputData, __inout InputTransformation &transform );

	void ConvertOutputData( __in size_t count, __in_ecount( count ) MeasurementSample *pSamples, __out TDFast *pOut );

} }