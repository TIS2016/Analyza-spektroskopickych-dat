#include "stdafx.h"
#include "StringUtils.h"

namespace DataAnalysis { namespace Utils {

	void LStrToStr ( __in LStrPtr pIn, __out std::string &outStr ) {
		unsigned char *pSrc = pIn->str;
		size_t charCount = static_cast<size_t>( pIn->cnt );
		for ( size_t i = 0; i <= charCount; i++ ) {
			outStr.push_back ( static_cast<char>( *pSrc ) );
			pSrc++;
		}
	}

	void RTrim ( __inout std::string &str, __in const char *pToTrim ) {
		size_t delPos = str.find_last_not_of ( pToTrim ) + 1;
		str.erase ( delPos );
	}

	void LTrim ( __inout std::string &str, __in const char *pToTrim ) {
		size_t delPos = str.find_first_not_of ( pToTrim );
		str.erase ( str.begin (), str.begin () + delPos );
	}

} }