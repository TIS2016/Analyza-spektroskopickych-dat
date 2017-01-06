#pragma once

#include "stdafx.h"
#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double, class SampleClass = MeasurementSample> class SummaryOperation : public IFunction<SampleClass> {
	public:
		void Initialize(
			__in const shared_ptr< IFunction<BaseType> > spPoly,
			__in const shared_ptr< IFunction<BaseType> > spTrig,
			__in const shared_ptr< IFunction<BaseType> > spSpline ) {

			if ( spPoly != nullptr && spTrig != nullptr && spSpline != nullptr ) {
				mSpPolynomial = spPoly;
				mSpTrigonometric = spTrig;
				mSpSpline = spSpline;
				mInitialized = true;
			}
		}

	protected:

		shared_ptr<IFunction<BaseType>> mSpPolynomial;
		shared_ptr<IFunction<BaseType>> mSpTrigonometric;
		shared_ptr<IFunction<BaseType>> mSpSpline;

	protected:

		BaseType GetSum( __in const BaseType &in ) const {
			BaseType sA, sB, sC;

			mSpPolynomial->Apply( in, sA );
			mSpTrigonometric->Apply( in, sB );
			mSpSpline->Apply( in, sC );

			return sA + sB + sC;
		}

	};

	template <class BaseType = double, class SampleClass = MeasurementSample> class SummaryDivision : public SummaryOperation<BaseType, SampleClass> {
	public:
		SummaryDivision() {
			mType = FT_SUMOP_DIV;
		};

		virtual inline void Apply( __in const SampleClass &in, __out SampleClass &out ) const {
			out.Y = in.Y / GetSum( in.X );
		}
	};

	template <class BaseType = double, class SampleClass = MeasurementSample> class SummaryMultiplication : public SummaryOperation<BaseType, SampleClass> {
	public:
		SummaryMultiplication() {
			mType = FT_SUMOP_MUL;
		};

		virtual inline void Apply( __in const SampleClass &in, __out SampleClass &out ) const {
			out.Y = in.Y * GetSum( in.X );
		}
	};

} }