#pragma once

#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

#pragma push_macro( "SPL_LIMITS" )
#ifdef min
	#undef min
#endif

#ifdef max
	#undef max
#endif

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
			Interpolates the spline at given point
			Should be done only between given control points
		*/
		virtual inline void Apply( __in const BaseType &in, __out BaseType &out ) const {
			Segment *pSegment = GetSegment( in );

			BaseType xScale = pSegment->pCpEnd->x - pSegment->pCpStart->x;
			BaseType inScaled = ( in - pSegment->pCpStart->x ) / xScale;
			BaseType S[4] = { inScaled*inScaled*inScaled, inScaled*inScaled, inScaled, 1 };
			BaseType C[4] = {   pSegment->pCpStart->fx,
								pSegment->pCpEnd->fx,
								pSegment->pCpStart->t * xScale,
								pSegment->pCpEnd->t * xScale };

			// P = S * h * C
			BaseType tmp[4] = { 0, 0, 0, 0 };
			Vect4MulMat44( S, m_basisFunctionsT, tmp );
			Vect4DotProduct( tmp, C, &out );
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
		
		const BaseType m_basisFunctionsT[16] = {  2, -2,  1,  1,
												 -3,  3, -2, -1,
												  0,  0,  1,  0,
												  1,  0,  0,  0 };

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

		BaseType GetLinearPoint( __in ControlPoint &pt0, __in ControlPoint &pt1, __in const BaseType xIn ) {
			BaseType k = ( pt1.fx - pt0.fx ) / ( pt1.x - pt0.x );
			BaseType q = pt1.fx - k * pt1.x;
			return k*xIn + q;
		}

		/*
			Convert input data to interal control points
			Input data format: [ CP1_X | CP1_Fx | CP2_X | CP2_Fx | .... | CPn_X | CPn_Fx ]
			
			Caller is responsible for:
				1) providing at least 2 control points to the function
				2) providing control points pre-sorted along X axis
				3) not including 2 control points with the same X value ( will cause problems in differentiation approximation )
		*/
		void ConvertToControlPoints( __in const size_t ptCount, __in_ecount( ptCount * 2 ) const BaseType *pVals ) {
			ControlPoint *pDst = m_controlPoints.Ptr() + 1;
			const BaseType *pSrc = pVals;
		
			for ( size_t i = 0; i < ptCount; i++ ) {
				pDst->x = *pSrc;
				pDst->fx = *(pSrc + 1);
				pSrc += 2;
				pDst++;
			}

			// "inf" and "-inf" points are actually first/last CP -/+ 100 times the span of the spline on X axis
			// they lie on a line connecting first and second CP ( or last and it's predecessor )
			BaseType xRange = ( *( pVals + 2 * ( ptCount - 1 ) ) ) - ( *pVals );
			m_controlPoints[0].x = ( *pVals - xRange * 100 );
			m_controlPoints[0].fx = GetLinearPoint( m_controlPoints[1], m_controlPoints[2], m_controlPoints[0].x );

			m_controlPoints[m_cpCount - 1].x = ( *( pVals + 2 * ( ptCount - 1 ) ) ) + xRange * 100;
			m_controlPoints[m_cpCount - 1].fx = GetLinearPoint( m_controlPoints[m_cpCount - 3], m_controlPoints[m_cpCount - 2], m_controlPoints[m_cpCount - 1].x );
		}


		void InitializeSegments() {
			ControlPoint *pCpSrc = m_controlPoints.Ptr();
			Segment *pDst = m_segments.Ptr();
			for ( size_t i = 0; i < m_segmentCount; i++ ) {
				pDst->pCpStart = pCpSrc++;
				pDst->pCpEnd = pCpSrc;
				pDst++;
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


		Segment* GetSegment( __in const BaseType &x ) const {
			// do binary search for correct segment
			size_t begin = 0;
			size_t end = m_segmentCount;
			
			while ( begin < end ) {
				size_t mid = begin + ( ( end - begin ) >> 1 );

				Segment *pSegment = m_segments.Ptr() + mid;
				if ( pSegment->pCpStart->x <= x )
				{
					if ( pSegment->pCpEnd->x > x ) {
						return pSegment;
					}
					else {
						begin = mid + 1;
					}
				}
				else {
					end = mid;
				}
 			}

			return ( m_segments.Ptr() + end );
		}
	};

#pragma pop_macro( "SPL_LIMITS" )

} }
