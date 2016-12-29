#pragma once

#include "stdafx.h"
#include "TransformationsLibPrivate.h"
#include "Polynomial.h"
#include "HermitePolynomialProvider.h"
#include "LegendrePolynomialProvider.h"

namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double> class IPolynomialTransform : public IFunction<BaseType> {
	public:
		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) BaseType *pCnsts ) = 0;

		void Initialize( __in const Buffer<BaseType> &params ) {
			Initialize( params.Length(), params.Ptr() );
		}

		void Initialize( __in const vector<BaseType> &params ) {
			Initialize( params.size(), params.data() );
		}
	};

	template <class BaseType = double> class BasicPolynomialTransform : public IPolynomialTransform<BaseType> {
	public:
		BasicPolynomialTransform() { 
			mType = FT_POLY_BASIC;
		};

		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) BaseType *pCnsts ) {
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

		virtual void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = BaseType( 0 );
			for ( uint exp = 0; exp <= mDegree; exp++ ) {
				out += mConstants[exp] * mPolynomials[exp]->GetFor( in );
			}
		}

	protected:

		Buffer< shared_ptr<Polynomial<BaseType>> > mPolynomials;
		Buffer< BaseType > mConstants;
		uint mDegree;

	protected:

		void InitializeInternal( __in const uint cnstCnt ) {
			mDegree = cnstCnt - 1;
			mConstants.Allocate( cnstCnt );
			mPolynomials.Allocate( cnstCnt );
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
		
		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) BaseType *pCnsts ) {
			InitializeInternal( cnstCount );
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
			mType = FT_POLY_CHEBYSHEV;
		}

		virtual void Initialize( __in const uint cnstCount, __in_ecount( cnstCount ) BaseType *pCnsts ) {
			InitializeInternal( cnstCount );
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