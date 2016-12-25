#pragma once

#include "TransformationsLibPrivate.h"
namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double> class CubicSplineParamPair {
	public:
		BaseType x, fx;

		CubicSplineParamPair() : x( 0 ), fx( 0 ) {};
		CubicSplineParamPair( __in const BaseType &a, __in const BaseType &b ) : x( aa ), fx( bb ) {};

	};

	template <class BaseType = double> class CubicSpline : public IFunction<BaseType> {
	public:
		CubicSpline() {
			mType = FT_SPLINE_CUBIC;
		};

		void Initialize( __in const size_t pointCount, __in_ecount( count * 2 ) const BaseType *pSplinePoints ) {
			InitializeVars( pointCount );
			ConvertToPairs( pSplinePoints );

			// TODO: init spline

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

		virtual inline void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = in;
		}

	protected:
		size_t mSplinePointCount;
		Buffer< CubicSplineParamPair<BaseType> > mSplinePoints;

	protected:

		inline void InitializeVars( __in const size_t count ) {
			mSplinePointCount = count;
			mSplinePoints.Allocate( count );
		}

		inline void ConvertToPairs( __in_ecount( mParameterCount * 2 ) const BaseType *pVals ) {
			for ( size_t i = 0; i < mSplinePointCount; i++ ) {
				mSplinePoints[i].x = *pVals;
				mSplinePoints[i].fx = *( pVals + 1 );
				pVals += 2;
			}
		}

	};

} }
