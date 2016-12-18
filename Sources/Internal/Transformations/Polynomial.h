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

		void Initialize( __in const size_t count, __in_ecount( count ) const BaseType *pParams ) {
			if ( count > 0 ) {
				mDegree = count - 1;
				mConstants.Set( count, pParams );
				mInitialized = true;
			}
		}

		void Initialize( __in const Buffer<BaseType> &params ) {
			Initialize( params.Length(), params.Ptr() );
		}

		void Initialize( __in const vector<BaseType> &params ) {
			Initialize( params.size(), params.data() );
		}

		/** 
		Implementation for basic polynomial function in form a*n^0 + b*n^1 + c*n^2 + ... + {alpha}*n^mDegree
		Using: 
			- mDegree as max degree of polynomial
			- mValues as {a, b, c, ... } parameters in various polynomial degrees
		*/
		virtual void Apply ( __in const BaseType &in, __out BaseType &out ) const {
			out = mConstants[0];
			for (size_t exp = 1; exp <= mDegree; exp++) {
				out += pow(in, exp) * mConstants[exp];
			} 
		};

	protected:
		UINT mDegree;
		Buffer<BaseType> mConstants;
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