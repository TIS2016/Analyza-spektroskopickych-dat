#include "stdafx.h"
#include "StringUtils.h"

BEGIN_NAMESPACE

void LStrToStr(__in LStrPtr pIn, __out std::string &outStr) {
	unsigned char *pSrc = pIn->str;
	size_t charCount = static_cast<size_t>(pIn->cnt);
	for (size_t i = 0; i <= charCount; i++) {
		outStr.push_back(static_cast<char>(*pSrc));
		pSrc++;
	}
}

END_NAMESPACE