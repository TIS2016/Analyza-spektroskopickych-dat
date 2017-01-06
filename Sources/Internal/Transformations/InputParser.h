#pragma once

#include "TransformationsLib.h"

namespace DataAnalysis { namespace Transformations {

	void ConvertToInternal(__in const TransformationHeader &in, __out TransformationHeaderInternal &out);

	shared_ptr<IFunction<MeasurementSample>> GetTransformation( __in const TransformationHeader &transformInfo, __in const InputTransformation *pCaller );

} }