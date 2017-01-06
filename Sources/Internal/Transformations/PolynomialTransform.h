#pragma once

#include "stdafx.h"
#include "TransformationsLibPrivate.h"
#include "Polynomial.h"
#include "HermitePolynomialProvider.h"
#include "LegendrePolynomialProvider.h"

namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double> class IPolynomialTransform : public IFunction<BaseType> {
		// dummy class for now, maybe one day, something will appear here
	};

	template <class BaseType = double> class BasicPolynomialTransform : public IPolynomialTransform<BaseType> {
	public:
		BasicPolynomialTransform() { 
			mType = FT_POLY_BASIC;
		};

		void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) const BaseType *pCnsts ) {
			mSpPolynomial = shared_ptr<Polynomial<BaseType>>( new Polynomial<BaseType>( cnstCount - 1, pCnsts ) );
			if ( mSpPolynomial != nullptr ) {
				mInitialized = true;
			}
		}

		/** 
		Implementation for basic polynomial function in form a*n^0 + b*n^1 + c*n^2 + ... + {alpha}*n^mDegree
		Using: 
			- mDegree as max degree of polynomial
			- mValues as {a, b, c, ... } parameters in various polynomial degrees
		*/
		virtual void Apply ( __in const BaseType &in, __out BaseType &out ) const {
			out = mSpPolynomial->GetFor( in );
		};

	protected:

		shared_ptr<Polynomial<BaseType>> mSpPolynomial;

	};


	/*
		Interface for composite polynomial transformations
		eg. P(x) = c0*P0(x) + c1*P1(x) + ... + cn*Pn(x)

		TODO:
			1) add rescaling of input argument of Apply function to interval [-1, 1] (somewhere, somehow)
	*/
	template <class BaseType = double> class ICompositePolynomialTransform : public IPolynomialTransform<BaseType> {

	public:

		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) const BaseType *pCnsts, __in const BaseType &sclMin, __in const BaseType &sclMax ) = 0;

		virtual void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = BaseType( 0 );
			BaseType inTransformed = GetValueLocal( in );
			for ( uint exp = 0; exp <= mDegree; exp++ ) {
				out += mConstants[exp] * mPolynomials[exp]->GetFor( inTransformed );
			}
			// ?? rescale value back or not ?? who knows ??
		}

	protected:

		Buffer< shared_ptr<Polynomial<BaseType>> > mPolynomials;
		Buffer< BaseType > mConstants;
		uint mDegree;

		BaseType scaleMin;
		BaseType scaleMax;
		BaseType scaleInv;
		BaseType scale;

	protected:

		void InitializeInternal( __in const uint cnstCnt ) {
			mDegree = cnstCnt - 1;
			mConstants.Allocate( cnstCnt );
			mPolynomials.Allocate( cnstCnt );
		}

		void InitializeScale( __in const BaseType &sclMin, __in const BaseType &sclMax ) {
			scaleMin = sclMin;
			scaleMax = sclMax;
			scale = sclMax - sclMin;
			scaleInv = 1 / scale;
		}

		/*
			Rescales given value to fit into interval [-1, 1]
		*/
		BaseType GetValueLocal( __in const BaseType &value ) const {
			return ( ( ( value - scaleMin ) * 2 ) * scaleInv ) - 1;
		}

		/*
			Rascales value back from [-1, 1]
		*/
		BaseType GetValueGlobal( __in const BaseType &value ) const {
			return ( ( ( value + 1 ) / 2 ) * scale ) + scaleMin;
		}

	};

	/*
	Composite polynomial tranformation, based on Legendre polynomials
	*/
	template <class BaseType = double> class LegendrePolynomialTransform : public ICompositePolynomialTransform<BaseType> {
	public:

		LegendrePolynomialTransform() {
			mType = FT_POLY_LEGENDRE;
		}
		
		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) const BaseType *pCnsts, __in const BaseType &sclMin, __in const BaseType &sclMax ) {
			InitializeInternal( cnstCount );
			InitializeScale( sclMin, sclMax );
			LegendrePolynomialProvider<BaseType> provider;

			for ( uint i = 0; i < cnstCount; i++ ) {
				mPolynomials[i] = provider.GetPolynomial( i );
				mConstants[i] = *pCnsts;
				pCnsts++;
			}

			mInitialized = true;
		}

	};

	/*
		Composite polynomial tranformation, based on (Physicist's) Hermite polynomials
	*/
	template <class BaseType = double> class HermitePolynomialTransform : public ICompositePolynomialTransform<BaseType> {
	public:

		HermitePolynomialTransform() {
			mType = FT_POLY_HERMITE;
		}

		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) const BaseType *pCnsts, __in const BaseType &sclMin, __in const BaseType &sclMax ) {
			InitializeInternal( cnstCount );
			InitializeScale( sclMin, sclMax );
			HermitePolynomialProvider<BaseType> provider;

			for ( uint i = 0; i < cnstCount; i++ ) {
				mPolynomials[i] = provider.GetPolynomial( i );
				mConstants[i] = *pCnsts;
				pCnsts++;
			}

			mInitialized = true;
		}

	};

} }