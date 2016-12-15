#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {
	
	/*
		Trignonometric function which looks like this:
			Out = a * trigF( 2PI *b * In + c ),
			where tigF is sin or cos
	*/
	template <class BaseType> class TrigonometricFunction : public IFunction<BaseType> {
	public:
		void Initialize( __in const BaseType a, __in const BaseType b, __in const BaseType c ) {
			mA = a;
			mB = b;
			mC = c;
			mInitialized = true;
		}

		void Initialize( __in const Buffer<BaseType> &params ) {
			if ( params.Length() >= 3 ) {
				mA = params[0];
				mB = params[1];
				mC = params[2];
				mInitialized = true;
			}
		}

		void Initialize( __in_ecount( 3 ) const BaseType *pParams ) {
			mA = *pParams;
			mB = *( pParams + 1 );
			mC = *( pParams + 2 );
			mInitialized = true;
		}

		void Initialize( __in const vector<BaseType> &params ) {
			if ( params.size() >= 3 ) {
				mA = params[0];
				mB = params[1];
				mC = params[2];
				mInitialized = true;
			}
		}

		/* TODO: redo
		static shared_ptr<IFunction> GetFunction( __in TRIGONOMETRIC_TYPE type, __in void *pVals ) {
			switch ( type ) {
			case TT_SIN:
				return shared_ptr<IFunction>( new SinFunction<BaseType>( *static_cast<Buffer<BaseType> *>( pVals ) ) );
			case TT_COS:
				return shared_ptr<IFunction>( new CosFunction<BaseType>( *static_cast<Buffer<BaseType> *>( pVals ) ) );
			default:
				return nullptr;
			}
		}
		*/

	protected:
		BaseType mA, mB, mC;
	};

	// Out = a * sin( In * 2PI*b + c )
	template <class BaseType = double > class SinFunction : public TrigonometricFunction<BaseType> {
	public:
		SinFunction() { mType = FT_TRIG_SIN; };

		virtual inline void Apply( __in const BaseType &in, __out BaseType &out  ) const {
			out = mA * sin( 2 * M_PI*mB*in + mC );
		};
	};

	// Out = a * cos( In * 2PI*b + c )
	template <class BaseType = double > class CosFunction : public TrigonometricFunction<BaseType> {
	public:
		CosFunction() { mType = FT_TRIG_COS; };

		virtual inline void Apply( __in const BaseType &in, __out BaseType &out ) const {
			out = mA * cos( 2 * M_PI*mB*in + mC );
		};
	};

} }