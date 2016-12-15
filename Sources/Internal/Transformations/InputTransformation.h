#pragma once

#include "TransformationsLib.h"

namespace DataAnalysis { namespace Transformations {

	/*	TODO:
			[unneccessary bullshit]
			1) think of a better name than FunctionHeader
			2) think of moving FunctionHeader somewhere else ( where it makes more sense, Lib maybe? )
			3) think of better names for FunctionHeader members
			[important stuff]
			4) adapt ImputTransformation to new template-based Functions
			5) add input parsing ( eg. FunctionHeader -> IFunction(s) )
				- right now, there are a few hard-coded methods, generic methods/class would be nice
	*/

	struct FunctionHeader {
		string name;
		Buffer<string> subFunctions;
		Buffer<Buffer<double>> functValues;
	};

	class InputTransformation {
	public:
		InputTransformation() {};

		void AddSample( __in MeasurementSample &sample );
		void AddSample( __in const double X, __in const double Y, __in const double uncertainty );

		void AddFunction( __in shared_ptr< IFunction<MeasurementSample> > spFunct );
		void AddFunction( __in FunctionHeader &funct ); 

		void CalculateTransformations( __out Buffer<MeasurementSample> &output );
		void CalculateTransformations( __in size_t count, __out_ecount( count ) MeasurementSample *pOutput );

	protected:

		shared_ptr< IFunction<MeasurementSample> > CreateXTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames );

		shared_ptr< IFunction<MeasurementSample> > CreateYTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames );

	protected:

		vector< MeasurementSample > mInputSamples;
		vector< shared_ptr< IFunction<MeasurementSample> > > mTranformations;

	};
} }

