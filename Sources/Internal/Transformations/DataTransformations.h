#pragma once

#include "stdafx.h"
#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

	// TODO: redo to empty constructor + Initialize methods template

	/*
		Transformation function for data on X axis
		Xout = Scale(Xin) + Offset
			- Scale: can be any function, on basic type, which implements IFunction interface
			- Offset: any number, shifts X axis by constant
	*/
	template <class BaseType = double> class XTransform : public IFunction<MeasurementSample> {
	public:

		XTransform () {
			mType = FT_TRANSFORM_X;
		};
		
		void Initialize( __in const BaseType off, __in const shared_ptr< IFunction<BaseType> > spScaleFunction ) {
			if ( spScaleFunction != nullptr ) {
				mOffset = off;
				mSpScaleFunction = spS :DcaleFunction;
				mInitialized = true;
			}
		}

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			mSpScaleFunction->Apply ( in.X, out.X );
			out.X += mOffset;
		}

	protected:
		BaseType mOffset;
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
			double poly, trig = 0;
			mSpSpline->Apply ( in.X, out.X );
			mSpPolynomial->Apply ( in.X, poly );
			mSpTrigonometric->Apply ( in.X, trig );
			out.X += poly + trig;
		}

	protected:
		shared_ptr< IFunction<BaseType> > mSpPolynomial;
		shared_ptr< IFunction<BaseType> > mSpTrigonometric;
		shared_ptr< IFunction<BaseType> > mSpSpline;
	};

} }