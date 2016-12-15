#pragma once

#include "TransformationsLibPrivate.h"
#include "Polynomial.h"
#include "Spline.h"
#include "Trigonometric.h"
#include "CompositeFunctions.h"
#include "DataTransformations.h"
/*
#include "InputTransformation.h"
*/

namespace DataAnalysis { namespace Transformations {

	// Function for converting internal function type enum values to input/output function types
	inline int GetInputFunctionType( __in const FUNCTION_TYPE type ) {
		int resType = static_cast<int>(type);
		while ( resType >= 0x1000 ) {
			resType -= 0x1000;
		}
		return resType;
	}

	/* 
		Function for retrieving internal function type based on information from labview input data
		Input:	functName: string indetifying sub-function
				typeId:	int indetifying concrete sub-function type
		Mapping:
			Polynomial functions: XScl, YOff, YPol
				Types:	0 -> basic polynomial
						1 -> Legendre polynomial
						2 -> Chebyshev polynomial
			Trigonometric functions: YTrg
				Types:	0 -> Sin
						1 -> Cos
			Spline functions: YSpl
				Types: None
			Summary operation: YTyp
				Types:	1 -> Multiplication
						2 -> Division
	*/
	inline FUNCTION_TYPE GetInternalFunctionType( __in const string &functName, __in int typeId = 0 ) {
		if ( functName.compare( "XScl" ) == 0 || functName.compare( "YOff" ) == 0 || functName.compare( "YPol" ) == 0 ) {
			typeId += 0x1000;
		}
		else if ( functName.compare( "YTrg" ) == 0 ) {
			typeId += 0x2000;
		}
		else if ( functName.compare( "YSpl" ) == 0 ) {
			typeId += 0x3000;
		}
		else if ( functName.compare( "YTyp" ) == 0 ) {
			typeId += 0x4000;
		}
		else {
			return FT_UNKNOWN;
		}

		return static_cast<FUNCTION_TYPE>( typeId );
	}

} }