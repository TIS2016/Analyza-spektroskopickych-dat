#pragma once

#include "stdafx.h"
#include "Polynomial.h"

namespace DataAnalysis { namespace Transformations {

	template <class PolyBaseType = double > class IRecursivePolynomialProvider
	{
	public:

		virtual ~IRecursivePolynomialProvider() {};

		virtual shared_ptr<Polynomial<PolyBaseType>> GetPolynomial( __in const uint degree ) = 0;

		friend class Polynomial<PolyBaseType>;

	protected:

		vector< shared_ptr< Polynomial<PolyBaseType> > > m_polynomials;

	};

} }