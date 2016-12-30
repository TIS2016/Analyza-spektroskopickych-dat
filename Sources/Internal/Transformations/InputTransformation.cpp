#include "stdafx.h"
#include "InputTransformation.h"
#include "InputParser.h"

namespace DataAnalysis { namespace Transformations {

	void InputTransformation::AddSample( __in const MeasurementSample &sample ) {
		mInputSamples.push_back( sample );
		UnsetFlag( ITS_NO_SAMPLES );
	};

	void InputTransformation::AddSample( __in const double x, __in const double y, __in const double uncertainty ) {
		mInputSamples.push_back( MeasurementSample( x, y, uncertainty ) );
		UnsetFlag( ITS_NO_SAMPLES );
	}

	void InputTransformation::AddTransformation( __in const shared_ptr< IFunction<MeasurementSample> > spTransform ) {
		FUNCTION_TYPE type = spTransform->GetType();

		switch ( type )
		{
		case( FT_MODEL_BASELINE ):
		case( FT_MODEL_PEAKS ):
			mModelTransformation.push_back( spTransform );
			UnsetFlag( ITS_NO_MODEL );
			break;
		case( FT_TRANSFORM_X ):
		case( FT_TRANSFORM_Y ):
			mTranformations.push_back( spTransform );
			UnsetFlag( ITS_NO_TRANSFOMATIONS );
			break;
		default:
			break;
		}

	};

	void InputTransformation::AddTransformation( __in const TransformationHeader &transform ) {
		shared_ptr<IFunction<MeasurementSample>> spTransform = GetTransformation( transform );
		if ( spTransform ) {
			AddTransformation( spTransform );
		}
	}

	void InputTransformation::CalculateTransformations( __in const size_t count, __out_ecount( count ) MeasurementSample *pOutput ) {
		size_t inputSize = mInputSamples.size();
		if ( FlagIsSet( ITS_NO_SAMPLES | ITS_NO_TRANSFOMATIONS ) || count < inputSize ) {
			return;
		}

		for ( auto transformIt = mTranformations.begin(); transformIt != mTranformations.end(); transformIt++ ) {
			shared_ptr<IFunction<MeasurementSample>> spTransformation = *transformIt;
			spTransformation->ApplyOnData( inputSize, mInputSamples.data(), pOutput );
		}

		SetFlag( ITS_TRANSFORM_DONE );
	}

	void InputTransformation::CalculateTransformations( __out Buffer<MeasurementSample> &output ) {
		size_t inputSize = mInputSamples.size();
		if ( output.Length() < inputSize ) {
			output.Allocate( inputSize );
		}
		
		CalculateTransformations( inputSize, output.Ptr() );
		// copy output data to internal, so that model can be calculated
		mOutputSamples.assign( output.Ptr(), output.Ptr() + inputSize );
	}

	void InputTransformation::CalculateTransformations() {
		size_t sampleCount = mInputSamples.size();
		mOutputSamples.clear();
		mOutputSamples.reserve( sampleCount );
		CalculateTransformations( sampleCount, mOutputSamples.data() );
	}

	inline bool InputTransformation::FlagIsSet( __in const int flag ) const {
		return ( mState & flag ) != 0;
	}

	inline void InputTransformation::SetFlag( __in const int flag ) {
		mState |= flag;
	}

	inline void InputTransformation::UnsetFlag( __in const int flag ) {
		mState &= ~flag;
	}

} }