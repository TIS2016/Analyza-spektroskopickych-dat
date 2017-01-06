#pragma once

namespace DataAnalysis { namespace Transformations {
	
	// Type converting functions (LW<->Internal) depend heavily on enum's structure
	enum FUNCTION_TYPE {
		FT_POLY_BASIC = 0x1000,
		FT_POLY_LEGENDRE = 0x1001,
		FT_POLY_HERMITE = 0x1002,

		FT_TRIG_SIN = 0x2000,
		FT_TRIG_COS = 0x2001,

		FT_SPLINE_CUBIC_HERMITE = 0x3000,

		FT_SUMOP_UNKNOWN = 0x4000,
		FT_SUMOP_DIV = 0x4001,
		FT_SUMOP_MUL = 0x4002,

		FT_CONSTANT = 0x5000,

		FT_TRANSFORM_X = 0x9000,
		FT_TRANSFORM_Y = 0x9001,
		FT_MODEL = 0x9002,
		FT_MODEL_BASELINE = 0x9003,
		FT_MODEL_PEAKS = 0x9004,

		FT_UNKNOWN = 0xFFFFFFFF
	};

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

	template <class BaseType> class IFunction {
		/*	
			SUGGESTIONS:
				1) HRESULT return type for Apply methods
				2) Initialize template with pass-by-value arguments
				3) Apply returning instead of filling up a reference 
		*/
	public:

		virtual ~IFunction() {};

		template < typename SubClass, typename ... ArgTypes > void Initialize( ArgTypes& ... args ) {
			static_cast<SubClass *>(this)->Initialize( args... );
		} 

		/*
		Proposed variant of Initialize method, which can be more comfortable to use in some cases.
		Pass by value means that arguments have to be passed by pointer to avoid copy-ing

		template < typename SubClass, typename ... ArgTypes > void Initialize( ArgTypes ... args ) {
			static_cast<SubClass *>(this)->Initialize( args... );
		}
		*/

		// for the sake of performance, do not perform initialization checks in these
		virtual void Apply ( __in const BaseType &in, __out BaseType &out ) const = 0;

		/*
		Proposed additional Apply method, easier and cleaner to use in many cases.
		Although, to avoid copy-ing of data structures, this would have to be used only on simple types
		
		As not all classes inheriting this interface will have a use for this function (eg. will not implement their own version), it's not pure virtual (abstract)

		virtual BaseType Apply( __in const BaseType in ) const {
			return BaseType(0);
		}
		*/

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