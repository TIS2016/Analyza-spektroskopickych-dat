#pragma once

#include "stdafx.h"

namespace DataAnalysis { namespace Transformations {

	/*
		Non-compact Array-based Polynomial class
	*/
	template <class BaseType = double> class Polynomial {
	public:
		// -- Construction
		Polynomial() {
			m_coefCnt = 0;
		};

		Polynomial( __in const size_t deg, __in_ecount(deg + 1) const BaseType *pCoefs ) {
			m_degree = deg;
			m_coefCnt = deg + 1;
			m_coefs.reserve( m_coefCnt );
			m_coefs.assign( pCoefs, pCoefs + m_coefCnt );
		};

		Polynomial( initializer_list<BaseType> coefs ) {
			m_degree = coefs.size() - 1;
			m_coefCnt = coefs.size();
			m_coefs.assign( coefs );
		}

		Polynomial( __in const shared_ptr<Polynomial<BaseType>> spSource ) {
			m_degree = spSource->m_degree;
			m_coefCnt = spSource->m_coefCnt;
			m_coefs = spSource->m_coefs;
		};

		// -- Modification
		void SetCoef( __in const uint deg, __in const BaseType coef ) {
			if ( m_coefCnt <= deg ) {
				Grow( deg + 1 );
			}

			m_coefs[deg] = coef;
		};

		void Multiply( __in const BaseType cnst ) {
			for ( size_t i = 0; i < m_coefCnt; i++ ) {
				m_coefs[i] *= cnst;
			}
		}

		void Multiply( __in const BaseType cnst, __in const uint deg ) {
			Multiply( cnst );
			m_coefCnt += deg;
			m_degree += deg;
			m_coefs.insert( m_coefs.begin(), deg, BaseType(0) );
		}
		 
		void Subtract( __in const shared_ptr<Polynomial<BaseType>> spPoly ) {
			size_t subLen = min( m_coefCnt, spPoly->m_coefCnt );
			for ( size_t i = 0; i < subLen; i++ ) {
				m_coefs[i] -= spPoly->m_coefs[i];
			}
		}

		// -- Usage

		uint GetDegree() const {
			return m_degree;
		}

		BaseType GetFor( __in const BaseType x ) const {
			BaseType res( 0 );
			for ( size_t exp = 0; exp < m_coefCnt; exp++ ) {
				res += m_coefs[exp] * static_cast<BaseType>( pow( x, exp ) );
			}
			return res;
		}


	protected:

		vector<BaseType> m_coefs;
		uint m_degree;
		size_t m_coefCnt;

	protected:

		void Grow( __in const size_t newCoefCnt ) {
			Buffer<BaseType> tmp;
			tmp.Set( m_coefs );

			m_coefs.resize( newCoefCnt, BaseType(0) );
			memcpy( m_coefs.data(), tmp.Ptr(), m_coefCnt * sizeof( BaseType ) );
			m_coefCnt = newCoefCnt;
		}

	};

} }