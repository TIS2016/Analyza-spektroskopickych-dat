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
		LegendrePolynomialProvider() {};
		
		~LegendrePolynomialProvider() {};

		virtual shared_ptr<Polynomial<BaseType>> GetPolynomial( __in const uint degree ) {
			UNREFERENCED_PARAMETER( degree );
			return nullptr;
		};

	protected:

	};

} }