#pragma once

#include "stdafx.h"
#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

	/*
		Transformation function for data on X axis
		Xout = Scale(Xin) + Offset(Xin)
			- Scale: can be any function, on basic type, which implements IFunction interface
			- Offset: any function on basic type implementing IFunction interface
	*/
	template <class BaseType = double> class XTransform : public IFunction<MeasurementSample> {
	public:

		XTransform () {
			mType = FT_TRANSFORM_X;
		};

		void Initialize( __in const shared_ptr<IFunction<BaseType>> spOffsetFunction, __in const shared_ptr< IFunction<BaseType> > spScaleFunction ) {
			if ( spOffsetFunction != nullptr && spScaleFunction != nullptr ) {
				mSpOffsetFunction = spOffsetFunction;
				mSpScaleFunction = spScaleFunction;
				mInitialized = true;
			}
		}

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			BaseType offset( 0 );
			mSpOffsetFunction->Apply( in.X, offset );
			mSpScaleFunction->Apply ( in.X, out.X );
			out.X += offset;
		}

	protected:

		shared_ptr< IFunction<BaseType> > mSpOffsetFunction; // will be constant function in our case
		shared_ptr< IFunction<BaseType> > mSpScaleFunction;

	};


	/*
		Transformation function for data on Y axis
		Yout = Offset(Xin) + Summary( SampleIn )
			- Offset: any function, on basic type, which implements IFunction interface ( polynomial in this case )
			- Summary: any funcion, on MeasurementSample type, which implements IFunction interface ( in our case: ( Yin *arithmetic op* ( Polynomial(Xin) + Trigonometric(Xin) + Spline(Xin) ) )
	*/
	template <class BaseType = double> class YTransform : public IFunction<MeasurementSample> {
	public:

		YTransform () {
			mType = FT_TRANSFORM_Y;
		};

		void Initialize( __in const shared_ptr< IFunction<BaseType> > spOffsetFunct, __in const shared_ptr< IFunction<MeasurementSample> > spSummaryFunct ) {
			if ( (spOffsetFunct != nullptr) && (spSummaryFunct != nullptr) ) {
				mSpOffsetFunction = spOffsetFunct;
				mSpSummaryFunction = spSummaryFunct;
				mInitialized = true;
			}
		}

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			BaseType YOff;
			mSpOffsetFunction->Apply ( in.X, YOff );
			mSpSummaryFunction->Apply ( in, out );
			out.Y += YOff;
		}

	protected:
		shared_ptr< IFunction<BaseType> > mSpOffsetFunction;
		shared_ptr< IFunction<MeasurementSample> > mSpSummaryFunction;
	};

	/*
		Function for Baseline calculation
		Baseline = Polynomial(Xout) + Trigonometric(Xout) + Spline(Xout)
			- Polynomial, Trigonometric, Spline: functions which implement IFunction interface on double

		Input: Samples on which XTransform has been already applied ( eg. XOut )
	*/
	template <class BaseType = double> class BaselineTransform : public IFunction<MeasurementSample> {
	public:
		BaselineTransform() {
			mType = FT_MODEL_BASELINE;
		};

		void Initialize( 
			__in const shared_ptr< IFunction<BaseType> > spPoly, 
			__in const shared_ptr< IFunction<BaseType> > spTrig, 
			__in const shared_ptr< IFunction<BaseType> > spSpline ) 
		{
			if ( ( spPoly != nullptr ) && ( spTrig != nullptr ) && ( spSpline != nullptr ) ) {
				mSpPolynomial = spPoly;
				mSpTrigonometric = spTrig;
				mSpSpline = spSpline;
				mInitialized = true;
			}
		}

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			BaseType poly(0), trig(0);
			mSpSpline->Apply ( in.X, out.Model );
			mSpPolynomial->Apply ( in.X, poly );
			mSpTrigonometric->Apply ( in.X, trig );
			out.Model += poly + trig;
		}

	protected:
		shared_ptr< IFunction<BaseType> > mSpPolynomial;
		shared_ptr< IFunction<BaseType> > mSpTrigonometric;
		shared_ptr< IFunction<BaseType> > mSpSpline;
	};

} }