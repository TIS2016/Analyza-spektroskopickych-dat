// DataAnalysisDebug.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

int main()
{
	printf("Trololo Wololo...Testujeme...\n");

	
	Buffer< LStr > in;
	in.Allocate(5);
	LStrPtr pIn = in.Ptr();
	for (int i = 0; i < 5; i++) {
		pIn->cnt = 4 - i;
		*(pIn->str) = static_cast<char>('a' + i);
		pIn++;
	}

	std::string out;
	DataAnalysis::Utils::LStrToStr(in.Ptr(), out);

	printf("%s", out.data() );

	DebugLogger logger(LOG_PATH_FILIP);
	logger.Log("pica blesky");
	
	

	int barrier{ 0 };
	scanf_s("%d", &barrier);

    return 0;
}

