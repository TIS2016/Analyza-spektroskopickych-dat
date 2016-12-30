#include "stdafx.h"
#include "StringUtils.h"

namespace DataAnalysis { namespace Utils {

	void LStrToStr ( __in const LStrPtr pIn, __out std::string &outStr ) {
		const char *pSrc = reinterpret_cast<const char *>( pIn->str );
		outStr.assign( pSrc );
	}

} }