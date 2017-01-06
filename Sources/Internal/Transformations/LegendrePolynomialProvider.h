#pragma once

#include "IRecursivePolynomialProvider.h"

namespace DataAnalysis { namespace Transformations {

	/*
		Provider for Legendre polynomials

		Generates polynomials of given degree based on this recursive formula:
		P_0(x) = 1
		P_1(x) = x
		P_i(x) = (1/i) * ( (2(i-1) + 1)x * P_[i-1](x) -  (i-1) * P_[i-2](x) )
	*/
	template <class BaseType = double> class LegendrePolynomialProvider: public IRecursivePolynomialProvider<BaseType>
	{
	public:
		LegendrePolynomialProvider() {
			AddInitialPolynomials();
		};
		
		~LegendrePolynomialProvider() {};

		virtual shared_ptr<Polynomial<BaseType>> GetPolynomial( __in const uint degree ) {
			if ( degree < m_polynomials.size() ) {
				return m_polynomials[degree];
			}

			// get P_i-1 and P_i-2 
			shared_ptr<Polynomial<BaseType>> spP1 = GetPolynomial( degree - 1 );
			shared_ptr<Polynomial<BaseType>> spP2 = GetPolynomial( degree - 2 );

			// (i-1) * P_[i-2]
			shared_ptr<Polynomial<BaseType>> spP2Multiplied( new Polynomial<BaseType>( spP2 ) );
			spP2Multiplied->Multiply( static_cast<BaseType> ( degree - 1  ) );

			// ( 2( i - 1 ) + 1 )x * P_[i - 1]
			shared_ptr<Polynomial<BaseType>> spP1Multiplied( new Polynomial<BaseType>( spP1 ) );
			spP1Multiplied->Multiply( static_cast<BaseType> ( ((degree - 1) << 1) + 1 ), 1 );

			// (spP1Multiplied - spP2Multiplied)
			spP1Multiplied->Subtract( spP2Multiplied );
			shared_ptr<Polynomial<BaseType>> spNew( new Polynomial<BaseType>( spP1Multiplied ) );
			
			// (1/i) * (spP1Multiplied - spP2Multiplied)
			spNew->Multiply( static_cast<BaseType> ( 1 / degree ) );

			// memoize new polynomial
			m_polynomials.push_back( spNew );

			return spNew;
		};

	protected:

		/*
		Adds base cases P0 and P1
		*/
		void AddInitialPolynomials() {
			//P_0(x) = 1
			shared_ptr< Polynomial<BaseType> > spP0( new Polynomial<BaseType>( { 1 } ) );

			// P_1( x ) = x
			shared_ptr< Polynomial<BaseType> > spP1( new Polynomial<BaseType>( { 0, 1 } ) );

			m_polynomials.push_back( spP0 );
			m_polynomials.push_back( spP1 );
		};

	};

} }