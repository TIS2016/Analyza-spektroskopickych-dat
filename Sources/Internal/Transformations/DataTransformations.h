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
	class XTransform : public IFunction<MeasurementSample> {
	public:

		XTransform ( __in double const off, __in const shared_ptr< IFunction<double> > spScaleFunction ) : mOffset ( off ), mSpScaleFunction ( spScaleFunction ) {};

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			mSpScaleFunction->Apply ( in.X, out.X );
			out.X += mOffset;
		}

	protected:
		const double mOffset;
		const shared_ptr< IFunction<double> > mSpScaleFunction;
	};


	/*
		Transformation function for data on Y axis
		Yout = Offset(Xin) + Summary( SampleIn )
			- Offset: any function, on basic type, which implements IFunction interface ( polynomial in this case )
			- Summary: any funcion, on MeasurementSample type, which implements IFunction interface ( in our case: ( Yin *arithmetic op* ( Polynomial(Xin) + Trigonometric(Xin) + Spline(Xin) ) )
	*/
	class YTransform : public IFunction<MeasurementSample> {
	public:

		YTransform ( __in const shared_ptr< IFunction<double> > spOffsetFunct, __in const shared_ptr< IFunction<MeasurementSample> > spSummaryFunct ) :
			mSpSummaryFunction ( spSummaryFunct ), mSpOffsetFunction ( spOffsetFunct ) {};

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			double YOff = 0;
			mSpOffsetFunction->Apply ( in.X, YOff );
			mSpSummaryFunction->Apply ( in, out );
			out.Y += YOff;
		}

	protected:
		const shared_ptr< IFunction<double> > mSpOffsetFunction;
		const shared_ptr< IFunction<MeasurementSample> > mSpSummaryFunction;
	};

	/*
		Function for Baseline calculation
		Baseline = Polynomial(Xout) + Trigonometric(Xout) + Spline(Xout)
			- Polynomial, Trigonometric, Spline: functions which implement IFunction interface on double

		Input: Samples on which XTransform has been already applied ( eg. XOut )
	*/
	class BaselineTransform : public IFunction<MeasurementSample> {

		BaselineTransform ( __in const shared_ptr< IFunction<double> > spPoly, __in const shared_ptr< IFunction<double> > spTrig, __in const shared_ptr< IFunction<double> > spSpline ) :
			mSpSpline ( spSpline ), mSpTrigonometric ( spTrig ), mSpPolynomial ( spPoly ) {};

		virtual inline void Apply ( __in const MeasurementSample &in, __out MeasurementSample &out ) const {
			double poly, trig = 0;
			mSpSpline->Apply ( in.X, out.X );
			mSpPolynomial->Apply ( in.X, poly );
			mSpTrigonometric->Apply ( in.X, trig );
			out.X += poly + trig;
		}

	protected:
		const shared_ptr< IFunction<double> > mSpPolynomial;
		const shared_ptr< IFunction<double> > mSpTrigonometric;
		const shared_ptr< IFunction<double> > mSpSpline;
	};

} }