#include "stdafx.h"
#include "InputParser.h"

namespace DataAnalysis { namespace Transformations {
	

}
void Transformations::ConvertToInternal(TransformationHeader & in, TransformationHeaderInternal & out) {
	out.transformationType = GetInternalFunctionType(in.name, 0);
	out.functionParameters = in.functValues;

	out.subFunctions.Allocate(in.subFunctions.Length);
	for (size_t i = 0; i < in.subFunctions.Length; i++){
		out.subFunctions[i] = GetInternalFunctionType(in.subFunctions[i], 0);
	}

}


}