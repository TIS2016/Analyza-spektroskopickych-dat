#pragma once

#include "TransformationsLib.h"

namespace DataAnalysis { namespace Transformations {
	
	/* TODO:



	*/

	shared_ptr<IFunction<MeasurementSample>> GetTransformation( __in string &transformTag, __in Buffer<string> subFuncts, __in Buffer<Buffer<double>> subFunctsParams );

	shared_ptr<IFunction<MeasurementSample>> GetTransformation( __in TransformationHeader &transformation );

	void ConvertToInternal(__in TransformationHeader &in, __out TransformationHeaderInternal &out);

} }