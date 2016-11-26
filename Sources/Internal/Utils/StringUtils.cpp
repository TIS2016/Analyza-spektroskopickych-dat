#include "stdafx.h"
#include "StringUtils.h"

BEGIN_NAMESPACE

void LStrToStr(__in LStrPtr pIn, __out std::string &outStr) {
	while (pIn->cnt >= 0) {
		outStr.push_back(*pIn->str);
		pIn++;
	}
}

END_NAMESPACE