#pragma once

#include "TransformationsLibPrivate.h"

namespace DataAnalysis { namespace Transformations {

	template <class BaseType = double> class ConstantFunction : public IFunction<BaseType> {
	public:

		ConstantFunction() {
			mType = FT_CONSTANT;
		}

		void Initialize( __in const BaseType cnst ) {
			mConstant = cnst;
			mInitialized = true;
		}

		virtual void Apply( __in const BaseType &in, __out BaseType &out ) const {
			UNREFERENCED_PARAMETER( in );
			out = mConstant;
		}

	protected:

		BaseType mConstant;

	};

} }