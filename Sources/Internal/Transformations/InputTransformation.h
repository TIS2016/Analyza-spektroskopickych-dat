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

	enum INPUT_TRANSFORMATION_STATE {
		ITS_NO_SAMPLES = 0x1,
		ITS_NO_TRANSFOMATIONS = 0x2,
		ITS_CALCULATION_DONE = 0x4
	};

	class InputTransformation {
	public:
		InputTransformation() { 
			mState = ITS_NO_SAMPLES | ITS_NO_TRANSFOMATIONS; 
		};

		void AddSample( __in const MeasurementSample &sample );
		void AddSample( __in const double x, __in const double y, __in const double uncertainty );

		void AddTransformation( __in const shared_ptr< IFunction<MeasurementSample> > spFunct );
		void AddTransformation( __in const TransformationHeader &funct ); 

		void CalculateTransformations( __in const size_t count, __out_ecount( count ) MeasurementSample *pOutput );
		void CalculateTransformations( __out Buffer<MeasurementSample> &output );
		void CalculateTransformations();

	protected:

		/*
		shared_ptr< IFunction<MeasurementSample> > CreateXTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames );

		shared_ptr< IFunction<MeasurementSample> > CreateYTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames );
		*/

		inline bool FlagIsSet( __in const int flag ) const;

		inline void SetFlag( __in const int flag );

		inline void UnsetFlag( __in const int flag );

	protected:

		int mState;

		vector< MeasurementSample > mInputSamples;
		vector< MeasurementSample > mOutputSamples;

		vector< shared_ptr< IFunction<MeasurementSample> > > mTranformations;
		shared_ptr< IFunction<MeasurementSample> > mModelTransformation;

	};
} }

