#include "stdafx.h"
#include "StringUtils.h"

namespace DataAnalysis { namespace Utils {

	void LStrToStr ( __in LStrPtr pIn, __out std::string &outStr ) {
		const char *pSrc = reinterpret_cast<const char *>( pIn->str );
		outStr.assign( pSrc );

		/*
		unsigned char *pSrc = pIn->str;
		size_t charCount = static_cast<size_t>( pIn->cnt );
		for ( size_t i = 0; i <= charCount; i++ ) {
			if ( *pSrc == 0 ) {
				break;
			}

			outStr.push_back ( static_cast<char>( *pSrc ) );
			pSrc++;
		}
		*/
	}

} }