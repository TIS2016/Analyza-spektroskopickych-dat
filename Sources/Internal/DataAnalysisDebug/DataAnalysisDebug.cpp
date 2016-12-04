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
		*(pIn->str) = 'a' + i;
		pIn++;
	}

	std::string out;
	LStrToStr(in.Ptr(), out);

	printf("%s", out.data() );

	DebugLogger logger("F:\\School\\3.rocnik\\Winter\\TIS\\Log01.txt");
	logger.Log("pica blesky");
	

	int barrier{ 0 };
	scanf_s("%d", &barrier);

    return 0;
}

