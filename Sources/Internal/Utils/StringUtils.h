#pragma once

namespace DataAnalysis { namespace Utils {

	void LStrToStr ( __in LStrPtr pIn, __out std::string &outStr );

	inline void RTrim( __inout std::string &str, __in const char* pToTrim = " \t\n\v\f\r" ) {
		size_t delPos = str.find_last_not_of( pToTrim ) + 1;
		str.erase( delPos );
	}

	inline void LTrim( __inout std::string &str, __in const char *pToTrim = " \t\n\v\f\r" ) {
		size_t delPos = str.find_first_not_of( pToTrim );
		str.erase( str.begin(), str.begin() + delPos );
	}

	inline void Split( __in std::string &str, __in const char chr, __out std::vector<string> &out ) {
		std::string tmp;
		for ( auto it = str.begin(); it != str.end(); it++ ) {
			if ( *it == chr ) {
				out.push_back( tmp );
				tmp.clear();
			}
			else {
				tmp.push_back( *it );
			}
		}
		
		if ( tmp.size() > 0 ) {
			out.push_back( tmp );
		}
	}

} }