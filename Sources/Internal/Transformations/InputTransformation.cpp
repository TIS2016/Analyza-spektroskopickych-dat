#include "stdafx.h"
#include "InputTransformation.h"

namespace DataAnalysis { namespace Transformations {

	void InputTransformation::AddSample( __in MeasurementSample &sample ) {
		mInputSamples.push_back( sample );
	};

	void InputTransformation::AddSample( __in const double X, __in const double Y, __in const double uncertainty ) {
		mInputSamples.push_back( MeasurementSample { X,Y,uncertainty } );
	}

	void InputTransformation::AddFunction( __in shared_ptr< IFunction<MeasurementSample> > spFunct ) {
		mTranformations.push_back( spFunct ); 
	};

	void InputTransformation::AddFunction( __in FunctionHeader &funct ) {
		if ( funct.name.compare( "XT" ) == 0 ) {
			mTranformations.push_back( CreateXTransform( funct.functValues, funct.subFunctions ) );
		} else if ( funct.name.compare( "YT" ) == 0 ) {
			mTranformations.push_back( CreateYTransform( funct.functValues, funct.subFunctions ) );
		} else {

		}
	}

	void InputTransformation::CalculateTransformations( __out Buffer<MeasurementSample> &output ) {
		size_t inputSize = mInputSamples.size();
		if ( output.Length() < inputSize ) {
			output.Allocate( inputSize );
		}
		
		for ( size_t transformI = 0; transformI < mTranformations.size(); transformI++ ) {
			shared_ptr<IFunction<MeasurementSample>> transform = mTranformations[transformI];
			transform->Apply( inputSize, mInputSamples.data(), output.Ptr() );
		}
	}
	void InputTransformation::CalculateTransformations( __in size_t count, __out_ecount( count ) MeasurementSample *pOutput ) {
		size_t inputSize = mInputSamples.size();
		if ( count < inputSize ) {
			return;
		}
		UNREFERENCED_PARAMETER( pOutput );
		// TODO
	}

	shared_ptr< IFunction<MeasurementSample> > InputTransformation::CreateXTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames ) {
		/* Expected:
			subFunctNames = [XOff, XScl]
			paramValues = [[off], [scl1, scl2]]
		*/
		UNREFERENCED_PARAMETER( subFunctNames );

		shared_ptr<IFunction<double>> spScaleFunct( new BasicPolynomial<double>( paramValues[1] ) );
		return shared_ptr< IFunction<MeasurementSample> >( new XTransform( paramValues[0][0], spScaleFunct ) );
	}

	shared_ptr< IFunction<MeasurementSample> > InputTransformation::CreateYTransform( __in Buffer<Buffer<double>> &paramValues, __in Buffer<string> &subFunctNames ) {
		/* Expected:
			subFunctNames = [YOFf, YTyp, YPol, YTrg, YSpl]
			paramValues = [[pType, p0, p1, p2,.., pn], [type], [pType, p0, p1, p2,...,pn], [tType, p0, p1, p2], [p00, p01, p10, p11, p20, p21,..., pn0, pn1]] 
		*/
		UNREFERENCED_PARAMETER( subFunctNames );

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

} }