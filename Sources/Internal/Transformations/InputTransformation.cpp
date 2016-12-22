#include "stdafx.h"
#include "InputTransformation.h"

namespace DataAnalysis { namespace Transformations {

	void InputTransformation::AddSample( __in const MeasurementSample &sample ) {
		mInputSamples.push_back( sample );
		UnsetFlag( ITS_NO_SAMPLES );
	};

	void InputTransformation::AddSample( __in const double x, __in const double y, __in const double uncertainty ) {
		mInputSamples.push_back( MeasurementSample( x, y, uncertainty ) );
		UnsetFlag( ITS_NO_SAMPLES );
	}

	void InputTransformation::AddTransformation( __in const shared_ptr< IFunction<MeasurementSample> > spFunct ) {
		FUNCTION_TYPE type = spFunct->GetType();

		if ( type == FT_MODEL ) {
			mModelTransformation = spFunct;
			UnsetFlag( ITS_NO_TRANSFOMATIONS );
		} else if ( type == FT_TRANSFORM_X || type == FT_TRANSFORM_Y ) {
			mTranformations.push_back( spFunct );
			UnsetFlag( ITS_NO_TRANSFOMATIONS );
		}
	};

	void InputTransformation::AddTransformation( __in const TransformationHeader &funct ) {
		UNREFERENCED_PARAMETER( funct );
		// shared_ptr< IFunction<MeasurementSample> > spTransform = GetTransformation( funct );
		// AddTransformation( spTransform );
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
	
	/*
	shared_ptr< IFunction<MeasurementSample> > InputTransformation::CreateXTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames ) {
		/* Expected:
			subFunctNames = [XOff, XScl]
			paramValues = [[off], [scl1, scl2]]
		
		UNREFERENCED_PARAMETER( subFunctNames );
		UNREFERENCED_PARAMETER( paramValues );

		shared_ptr<IFunction<double>> spScaleFunct( new BasicPolynomial<double>( paramValues[1] ) );
		return shared_ptr< IFunction<MeasurementSample> >( new XTransform( paramValues[0][0], spScaleFunct ) );
	}
	
	shared_ptr< IFunction<MeasurementSample> > InputTransformation::CreateYTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames ) {
		/* Expected:
			subFunctNames = [YOFf, YTyp, YPol, YTrg, YSpl]
			paramValues = [[pType, p0, p1, p2,.., pn], [type], [pType, p0, p1, p2,...,pn], [tType, p0, p1, p2], [p00, p01, p10, p11, p20, p21,..., pn0, pn1]] 
		
		UNREFERENCED_PARAMETER( subFunctNames );
		UNREFERENCED_PARAMETER( paramValues );

		Buffer<double> offsetPolynomialValues( paramValues[0].Length() - 1, paramValues[0].Ptr() + 1 );
		shared_ptr< IFunction<double> > spOffset = Polynomial<double>::GetPolynomial( static_cast<POLYNOMIAL_TYPE>( (int)paramValues[0][0] ), &offsetPolynomialValues );

		Buffer<double> polyValues( paramValues[2].Length() - 1, paramValues[2].Ptr() + 1 );
		shared_ptr< IFunction<double> > spPoly = Polynomial<double>::GetPolynomial( static_cast<POLYNOMIAL_TYPE>( (int)paramValues[2][0] ), &polyValues );

		Buffer<double> trigValues( paramValues[3].Length() - 1, paramValues[3].Ptr() + 1 );
		shared_ptr< IFunction<double> > spTrig = TrigonometricFunction<double>::GetFunction( static_cast<TRIGONOMETRIC_TYPE>( (int)paramValues[3][0] ), &trigValues );

		Buffer<CubicSplineParamPair<double>> splineValues;
		size_t pairCount = paramValues[4].Length() / 2;
		splineValues.Allocate( pairCount );
		double *pSrc = paramValues[4].Ptr();
		for ( size_t i = 0; i < pairCount; i++ ) {
			splineValues[i].a = *pSrc;
			splineValues[i].b = *( pSrc + 1 );
			pSrc += 2;
		}

		shared_ptr< IFunction<double> > spSpline( new CubicSpline<double>( splineValues ) );

		shared_ptr< IFunction<MeasurementSample> > spSummary = GetSummaryOperationFunction( static_cast<SUMMARY_OPERATION_TYPE>( (int)paramValues[1][0] ), spPoly, spTrig, spSpline );
		
		return shared_ptr< IFunction<MeasurementSample> >( new YTransform( spOffset, spSummary ) );
	}
	*/

	inline bool InputTransformation::FlagIsSet( __in const int flag ) const {
		return (mState & flag) != 0;
	}

	inline void InputTransformation::SetFlag( __in const int flag ) {
		mState |= flag;
	}

	inline void InputTransformation::UnsetFlag( __in const int flag ) {
		mState &= ~flag;
	}

} }