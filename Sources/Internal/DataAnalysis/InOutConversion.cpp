#include "stdafx.h"
#include "InOutConversion.h"

namespace DataAnalysis { namespace InputOutput {

	inline void GetTransformationName( __in LStrHandle *pName, __out FunctionHeader *pHeader ) {
		LStrToStr( **pName, pHeader->name );
		LTrim( pHeader->name );
		RTrim( pHeader->name );
	}

	inline void SplitSubfunctions( __in string &subFunctString, __out FunctionHeader *pHeader ) {
		vector<string> subFunctions;
		Split( subFunctString, '@', subFunctions );

		pHeader->subFunctions.Allocate( subFunctions.size() );
		pHeader->functValues.Allocate( subFunctions.size() );
		memcpy( pHeader->subFunctions.Ptr(), subFunctions.data(), subFunctions.size() * sizeof( string ) ); // TODO: Add set method to Buffer
	}

	inline void GetTransformationSubfunctions( __in LStrHandle *pString, __out FunctionHeader *pHeader ) {
		string tmp;
		LStrToStr( **pString, tmp );
		LTrim( tmp );
		RTrim( tmp );

		SplitSubfunctions( tmp, pHeader );
	}

	void GetTransformationHeaders( __in size_t count, __in_ecount( count ) TD2 *pSrcStruct, __out_ecount( count ) FunctionHeader *pDst ) {
		LStrHandle *pSrc = pSrcStruct->String;

		for ( size_t i = 0; i < count; i++ ) {

			GetTransformationName( pSrc, pDst );
			GetTransformationSubfunctions( pSrc + 1, pDst );

			pSrc += 2;
			pDst++;
		}
	}

	void GetFunctionParameters( __in TD4 *pParameterIndexes, __in TD3 *pParameters, __inout Buffer<FunctionHeader> &headers ) {
		size_t transformParamAddressSpace = pParameterIndexes->dimSizes[1] * pParameterIndexes->dimSizes[2];
		size_t subFunctionParamAddressSpace = pParameterIndexes->dimSizes[2];

		size_t transformtionCount = headers.Length();

		int *pInedxes = pParameterIndexes->Numeric;

		for ( size_t transformI = 0; transformI < transformtionCount; transformI++ ) {
			FunctionHeader &transform = headers[transformI];
			size_t subFunctionCount = transform.subFunctions.Length();

			int *pIx = pInedxes;

			for ( size_t subFunctI = 0; subFunctI < subFunctionCount; subFunctI++ ) {

				size_t paramCount = *pIx;
				Buffer<double> &subfunctParams = transform.functValues[subFunctI];
				subfunctParams.Allocate( paramCount );
				for ( size_t paramIndexI = 0; paramIndexI < paramCount; paramIndexI++ ) {
					int index = *( pIx + paramIndexI + 1 );
					double value = *( pParameters->Numeric + 3 * index );
					subfunctParams[paramIndexI] = value;
				}

				pIx += subFunctionParamAddressSpace;
			}

			pInedxes += transformParamAddressSpace;
		}
	}

	void ConvertInputData( __in TD7 *pInputData, __inout InputTransformation &transform ) {
		size_t dataCount = ( *pInputData->X_in )->dimSize;

		double *pXin = ( *pInputData->X_in )->Numeric;
		double *pYin = ( *pInputData->Y_in )->Numeric;
		double *pDeviationIn = ( *pInputData->W_in )->Numeric;

		for ( size_t i = 0; i < dataCount; i++ ) {
			transform.AddSample( *pXin, *pYin, *pDeviationIn );

			pXin++;
			pYin++;
			pDeviationIn++;
		}
	}

	void ConvertOutputData( __in size_t count, __in_ecount( count ) MeasurementSample *pSamples, __out TDFast *pOut ) {
		double *pXout = ( *pOut->X_out )->Numeric;
		double *pYout = ( *pOut->Y_out )->Numeric;
		double *pDevOut = ( *pOut->W_out )->Numeric;
		MeasurementSample *pSrc = pSamples;

		for ( size_t sampleI = 0; sampleI < count; sampleI++ ) {
			*pXout = pSrc->X;
			*pYout = pSrc->Y;
			*pDevOut = pSrc->Deviation;

			pSrc++;
			pXout++;
			pYout++;
			pDevOut++;
		}
	}

} }
