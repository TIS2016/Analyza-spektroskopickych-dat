#include "stdafx.h"
#include "EntryPoint.h"
#include "InOutConversion.h"
#include <chrono>

using namespace DataAnalysis::InputOutput;

extern "C" _declspec(dllexport) int32_t fdata_fast(TD1 *DataPARin, TD7 *DataIN, TDFast *DataOUT_F, TD10 *Error) {
	UNREFERENCED_PARAMETER( Error );

	chrono::high_resolution_clock::time_point execStart = chrono::high_resolution_clock::now();

	TD1 *pInputParametes = DataPARin;
	TD7 *pInputData = DataIN;
	TDFast *pOutputData = DataOUT_F;
	
	InputTransformation transform;
	ConvertInputData( pInputData, transform );

	size_t transformationCount = ( *( pInputParametes->FuncNames ) )->dimSizes[0];

	Buffer<TransformationHeader> transformationHeaders;
	transformationHeaders.Allocate( transformationCount );
	GetTransformationHeaders( transformationCount, *pInputParametes->FuncNames, transformationHeaders.Ptr() );

	GetFunctionParameters( *pInputParametes->FuncParAdresses, *pInputParametes->ParamNumbers, transformationHeaders );
	
	for ( size_t i = 0; i < transformationHeaders.Length(); i++ ) {
		transform.AddTransformation( transformationHeaders[i] );
	}

	Buffer<MeasurementSample> output;
	transform.CalculateTransformations( output );
	transform.CalculateModel( output );

	ConvertOutputData( output.Length(), output.Ptr(), pOutputData );
	
	chrono::high_resolution_clock::time_point execEnd = chrono::high_resolution_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>( execEnd - execStart ).count();

	return static_cast<int32_t>( elapsed );
}

extern "C" _declspec(dllexport) int32_t fdata_complete(TD1 *DataPARin, TD7 *DataIN, TDComplete *DataOUT_C, TD10 *Error) {
	UNREFERENCED_PARAMETER(DataPARin);
	UNREFERENCED_PARAMETER(DataIN);
	UNREFERENCED_PARAMETER(DataOUT_C);
	UNREFERENCED_PARAMETER(Error);

	return 88;
}