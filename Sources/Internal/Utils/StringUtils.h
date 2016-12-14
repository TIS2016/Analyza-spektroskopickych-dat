#pragma once

namespace DataAnalysis { namespace Utils {

	void LStrToStr ( __in LStrPtr pIn, __out std::string &outStr );

	void RTrim ( __inout std::string &str, __in const char* toTrim );

	void LTrim ( __inout std::string &str, __in const char *pToTrim );

} }