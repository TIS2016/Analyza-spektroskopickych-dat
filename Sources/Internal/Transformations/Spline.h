#pragma once

#include "TransformationsLibPrivate.h"
namespace DataAnalysis { namespace Transformations {

#pragma push_macro( "SPL_INF" )
#ifdef INF
	#undef INF
#endif

#define INF 1.79769e+308

	template <class BaseType = double> class HermiteCubicSpline : public IFunction<BaseType> {
	public:
		HermiteCubicSpline() {
			mType = FT_SPLINE_CUBIC_HERMITE;
		}

		void Initialize( __in const size_t cpCount, __in_ecount( cpCount * 2 ) const BaseType *pControlPts ) {
			InitializeInternal( cpCount );
			ConvertToControlPoints( cpCount, pControlPts );
			InitializeSegments();
			InitializeDifferentiations();

			mInitialized = true;
		}

		void Initialize( __in const Buffer<BaseType> &params ) {
			// check if length is even
			if ( ( params.Length() & 1 ) == 0 ) {
				Initialize( params.Length() >> 1, params.Ptr() );
			}
		}

		void Initialize( __in const vector<BaseType> &params ) {
			// check if length is even
			if ( ( params.size() & 1 ) == 0 ) {
				Initialize( params.size() >> 1, params.data() );
			}
		}

		/*
			Sampels the spline at given point
		*/
		virtual inline void Apply( __in const BaseType &in, __out BaseType &out ) const {
			Segment *pSegment = GetSegment( in );
			BaseType xk = pSegment->pCpStart->x;
			BaseType xk1 = pSegment->pCpEnd->x;
			BaseType mk = pSegment->pCpStart->t;
			BaseType mk1 = pSegment->pCpEnd->t;

			BaseType t = ( in - xk ) / ( xk1 - in );
			BaseType vect[4] = { t, t, t, 1 };
			BaseType resVect[4] = { 0, 0, 0, 0 };
			Mat44MulVect4( m_basisFunctions, vect, resVect );

			out = ( resVect[0] * xk ) + ( resVect[1] * (xk1 - xk) * mk ) + ( resVect[2] * xk1 ) + ( resVect[3] * (xk1 - xk) * mk1 );
		}

	protected:

		struct ControlPoint {
			BaseType x;
			BaseType t; // tangent ( differentiation )
			BaseType fx;
		};

		struct Segment {
			ControlPoint *pCpStart;
			ControlPoint *pCpEnd;
		};


		size_t m_segmentCount;
		size_t m_cpCount;
		Buffer< ControlPoint > m_controlPoints;
		Buffer< Segment > m_segments;

		BaseType m_basisFunctions[16] = { 2, -3, 0, 1,
										  1, -2, 1, 0,
										  -2, 3, 0, 0,
										  1, -1, 0, 0 };

	protected:

		/*
		Allocates buffer for holding spline segment data
		There are ptCount + 1 segments, because we do include segments [-inf, CP1_X] and [CPn_X, inf]
		*/
		void InitializeInternal( __in const size_t ptCount ) {
			m_cpCount = ptCount + 2;
			m_controlPoints.Allocate( ptCount + 2 );
			m_segments.Allocate( ptCount + 1 );
			m_segmentCount = ptCount + 1;
		}

		/*
			Convert input data to interal control points
			Input data format: [ CP1_X | CP1_Fx | CP2_X | CP2_Fx | .... | CPn_X | CPn_Fx ]
			
			Caller is responsible for providing at least 2 control points to the function
		*/
		void ConvertToControlPoints( __in const size_t ptCount, __in_ecount( ptCount * 2 ) const BaseType *pVals ) {
			ControlPoint *pDst = m_controlPoints.Ptr();

			pDst->x = static_cast<BaseType>( -INF );
			if ( *( pVals + 1 ) < *( pVals + 3 ) ) {
				pDst->fx = static_cast<BaseType>( -INF );
			}
			else {
				pDst->fx = static_cast<BaseType>( INF );
			}
			pDst++;
			
		
			for ( size_t i = 0; i < ptCount; i++ ) {
				pDst->x = pVals;
				pDst->fx = pVals + 1;
				pVals += 2;
				pDst++;
			}

			pDst->x = static_cast<BaseType>( INF );
			if ( *( pVals - 3 ) < *( pVals - 1 ) ) {
				pDst->fx = static_cast<BaseType>( INF );
			}
			else {
				pDst->fx = static_cast<BaseType>( -INF );
			}
		}


		void InitializeSegments() {
			ControlPoint *pCpSrc = m_controlPoints.Ptr();
			Segment *pDst = m_segments.Ptr();
			for ( size_t i = 0; i < m_segmentCount; i++ ) {
				pDst->pCpStart = pCpSrc++;
				pDst->pCpEnd = pCpSrc;
			}
		}

		/*
			Basic differentiation approximation based on 2 neighbouring control points
						f(x_[i+1]) - f(x_[i-1])
			f'(x_i) = ----------------------------
							x_[i+1] - x_[i-1]
		*/
		inline BaseType ApproximateDifferentiation( ControlPoint *pPt0, ControlPoint *pPt2 ) {
			return ( pPt2->fx - pPt0->fx ) / ( pPt2->x - pPt0->x );
		}

		void InitializeDifferentiations() {
			ControlPoint *pPoint = m_controlPoints.Ptr() + 1;

			// skip "-inf" and "inf" CPs
			for ( size_t i = 0; i < m_cpCount - 2; i++ ) {
				pPoint->t = ApproximateDifferentiation( pPoint - 1, pPoint + 1 );
				pPoint++;
			}

			// "inf"/"-inf" CPs will have tanget same as CP0/CPn (linear)
			m_controlPoints[0].t = m_controlPoints[1].t;
			m_controlPoints[m_cpCount - 1].t = m_controlPoints[m_cpCount - 2].t;
		}


		Segment* GetSegment( __in const BaseType &x ) {
			// do binary search for correct segment
			size_t begin = 0;
			size_t end = m_segmentCount;
			
			while ( begin < end ) {
				size_t mid = ( end - begin ) >> 1;

				Segment *pSegment = m_segments.Ptr() + mid;
				if ( pSegment->pCpStart->x <= x )
				{
					if ( pSegment->pCpEnd->x >= x ) {
						return pSegment;
					}
					else {
						start = mid + 1;
					}
				}
				else {
					end = mid;
				}
 			}
		}
	};

#pragma pop_macro("SPL_INF")

} }
