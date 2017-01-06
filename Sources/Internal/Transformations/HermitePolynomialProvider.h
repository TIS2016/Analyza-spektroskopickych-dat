#pragma once

#include "stdafx.h"
#include "IRecursivePolynomialProvider.h"


namespace DataAnalysis { namespace Transformations {

	/*
		Provider for (Physicist's) Hermite polynomials
		
		Generates polynomials of given degree based on this recursive formula:
			H_0(x) = 1
			H_1(x) = 2x
			H_i(x) = 2x * H_[i-1](x) - 2(i-1) * H_[i-2](x)
	*/
	template< class BaseType = double > class HermitePolynomialProvider: public IRecursivePolynomialProvider<BaseType>
	{
	public:
		HermitePolynomialProvider() {
			AddInitialPolynomials();
		};

		~HermitePolynomialProvider() {};

		virtual shared_ptr<Polynomial<BaseType>> GetPolynomial( __in const uint degree ) {
			if ( degree < m_polynomials.size() ) {
				return m_polynomials[degree];
			}

			// get H_i-1 and H_i-2 
			shared_ptr<Polynomial<BaseType>> spH1 = GetPolynomial( degree - 1 );
			shared_ptr<Polynomial<BaseType>> spH2 = GetPolynomial( degree - 2 );
			
			// H_i = 2x * H_[i-1] - 2(i - 1) * H_[i-2]
			shared_ptr<Polynomial<BaseType>> spH2Multiplied( new Polynomial<BaseType>( spH2 ) );
			spH2Multiplied->Multiply( static_cast<BaseType> ( ( degree - 1 ) << 1 ) );

			shared_ptr<Polynomial<BaseType>> spNew( new Polynomial<BaseType>( spH1 ) );
			spNew->Multiply( 2, 1 );

			spNew->Subtract( spH2Multiplied );
			
			// memoize new polynomial
			m_polynomials.push_back( spNew );

			return spNew;
		};

	protected:

		/*
			Adds base cases H0 and H1
		*/
		void AddInitialPolynomials() {
			//H_0(x) = 1
			shared_ptr< Polynomial<BaseType> > spH0( new Polynomial<BaseType>( { 1 } ) );

			// H_1( x ) = 2x
			shared_ptr< Polynomial<BaseType> > spH1( new Polynomial<BaseType>( { 0, 2 } ) );

			m_polynomials.push_back( spH0 );
			m_polynomials.push_back( spH1 );
		};

	};

} }