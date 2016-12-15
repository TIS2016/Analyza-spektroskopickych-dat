#pragma once

#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

	// TODO: redo templates to empty constructor + Initialize methods ( like in Trigonometric.h, Spline.h )
	//			-> need to know what Legendre and Chebyshev polynomial look like ( their formula, parameters,.. )

	template <class BaseType = double> class Polynomial : public IFunction<BaseType> {	
	public:

		/* TODO: redo 
		static shared_ptr<IFunction> GetPolynomial( __in POLYNOMIAL_TYPE type, __in void *pVals ) {
			switch ( type ) {
			case PT_BASIC:
				return shared_ptr<IFunction>( new BasicPolynomial<BaseType>( *static_cast< Buffer<BaseType> * >( pVals ) ) );
				
				case: PT_LEGENDRE:
				return LegendrePolynomial<BaseType> ( type );
				case: PT_CHEBYSHEV:
				return ChebyshevPolynomial<BaseType> ( type );
				
			default:
				return nullptr;
			}
		}
		*/
	};


	template <class BaseType = double> class BasicPolynomial : public Polynomial<BaseType> {
	public:
		BasicPolynomial() { 
			mType = FT_POLY_BASIC;
		};

		virtual void Apply ( __in const BaseType &in, __out BaseType &out ) const {
			out = in;
		};

	protected:
		UINT mDegree;
		Buffer<BaseType> mValues;
	};

	template <class BaseType = double> class LegendrePolynomial : public Polynomial<BaseType> {
	public:

		LegendrePolynomial() {
			mType = FT_POLY_LEGENDRE;
		}

		virtual void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = in;
		};

	};

	template <class BaseType = double> class ChebyshevPolynomial : public Polynomial<BaseType> {
	public:

		ChebyshevPolynomial() {
			mType = FT_POLY_CHEBYSHEV;
		}

		virtual void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = in;
		};

	};

} }