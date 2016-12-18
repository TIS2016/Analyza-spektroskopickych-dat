#pragma once

namespace DataAnalysis { namespace Transformations {
	
	struct TransformationHeader {
		string name;
		Buffer<string> subFunctions;
		Buffer<Buffer<double>> functValues;
	};

	struct TransformationHeaderInternal {
		FUNCTION_TYPE transformationType;
		Buffer<FUNCTION_TYPE> subFunctions;
		Buffer<Buffer<double>> functionParameters;
	};

	// Type converting functions (LW<->Internal) depend heavily enum's structure
	enum FUNCTION_TYPE {
		FT_POLY_BASIC = 0x1000,
		FT_POLY_LEGENDRE = 0x1001,
		FT_POLY_CHEBYSHEV = 0x1002,

		FT_TRIG_SIN = 0x2000,
		FT_TRIG_COS = 0x2001,

		FT_SPLINE_CUBIC = 0x3000,

		FT_SUMOP_UNKNOWN = 0x4000,
		FT_SUMOP_DIV = 0x4001,
		FT_SUMOP_MUL = 0x4002,

		FT_TRANSFORM_X = 0x5000,
		FT_TRANSFORM_Y = 0x5001,
		FT_MODEL = 0x5002,
		FT_MODEL_BASELINE = 0x5003,
		FT_MODEL_PEAKS = 0x5004,

		FT_UNKNOWN = 0xFFFFFFFF
	};

	template <class BaseType> class IFunction {
		/*	TODO: 
				1) add function status enum
			
			SUGGESTIONS:
				1) HRESULT return type for Apply methods
		*/
	public:

		virtual ~IFunction() {};

		template < typename SubClass, typename ... ArgTypes > void Initialize( ArgTypes& ... args ) {
			static_cast<SubClass *>(this)->Initialize( args... );
		} 

		// for the sake of performance, do not perform initialization checks in these
		virtual void Apply ( __in const BaseType &in, __out BaseType &out ) const = 0;

		void ApplyOnData( __in const size_t count, __in_ecount( count ) const BaseType *pIn, __out_ecount( count ) BaseType *pOut ) const {
			if ( !mInitialized ) {
				return;
			}

			const BaseType *pSampleIn = pIn;
			BaseType *pSampleOut = pOut;
			for ( size_t sampleI = 0; sampleI < count; sampleI++ ) {
				Apply( *pSampleIn, *pSampleOut );
				pSampleIn++;
				pSampleOut++;
			}
		}

		FUNCTION_TYPE GetType() const { return mType; }

		bool GetStatus() const { return mInitialized; }

	protected:
		
		bool mInitialized = false;

		FUNCTION_TYPE mType;

	};

} }