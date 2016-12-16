#pragma once

#include "TransformationsLibPrivate.h"
namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double> class CubicSplineParamPair {
	public:
		BaseType a, b;

		CubicSplineParamPair() : a( 0 ), b( 0 ) {};
		CubicSplineParamPair( __in const BaseType &aa, __in const BaseType &bb ) : a( aa ), b( bb ) {};

	};

	template <class BaseType = double> class CubicSpline : public IFunction<BaseType> {
	public:
		CubicSpline() {
			mType = FT_SPLINE_CUBIC;
		};

		void Initialize( __in const Buffer<CubicSplineParamPair<BaseType>> &params ) {
			mParameterCount = params.Length();
			mParameters.Set( params );
			mInitialized = true;
		}

		void Initialize( __in const vector<CubicSplineParamPair<BaseType>> &params ) {
			mParameterCount = params.size();
			mParameters.Set( params );
			mInitialized = true;
		}

		void Initialize( __in const Buffer<BaseType> &params ) {
			// check if length is even
			if ( ( params.Length() & 1 ) == 0 ) {
				InitializeVars( params.Length() >> 1 );
				ConvertToPairs( params.Ptr() );

				mInitialized = true;
			}
		}

		void Initialize( __in const vector<BaseType> &params ) {
			// check if length is even
			if ( ( params.size() & 1 ) == 0 ) {
				InitializeVars( params.size() >> 1 );
				ConvertToPairs( params.data() );

				mInitialized = true;
			}
		}

		virtual inline void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = in;
		}

	protected:
		size_t mParameterCount;
		Buffer< CubicSplineParamPair<BaseType> > mParameters;

	protected:

		inline void InitializeVars( __in const size_t count ) {
			mParameterCount = count;
			mParameters.Allocate( count );
		}

		inline void ConvertToPairs( __in_ecount( mParameterCount * 2 ) const BaseType *pVals ) {
			for ( size_t i = 0; i < mParameterCount; i++ ) {
				mParameters[i].a = *pVals;
				mParameters[i].b = *( pVals + 1 );
				pVals += 2;
			}
		}

	};

} }
