#pragma once

#include "EntryPoint.h"

namespace DataAnalysis { namespace InputOutput {

	inline void GetTransformationName( __in LStrHandle *pName, __out TransformationHeader *pHeader );

	inline void SplitSubfunctions( __in string &subFunctString, __out TransformationHeader *pHeader );

	inline void GetTransformationSubfunctions( __in LStrHandle *pString, __out TransformationHeader *pHeader );

	void GetTransformationHeaders( __in size_t count, __in_ecount( count ) TD2 *pSrcStruct, __out_ecount( count ) TransformationHeader *pDst );

	void GetFunctionParameters( __in TD4 *pParameterIndexes, __in TD3 *pParameters, __inout Buffer<TransformationHeader> &headers );

	void ConvertInputData( __in TD7 *pInputData, __inout InputTransformation &transform );

	void ConvertOutputData( __in size_t count, __in_ecount( count ) MeasurementSample *pSamples, __out TDFast *pOut );

} }