// DataAnalysisDebug.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "../Transformations/TransformationsLib.h"
using namespace DataAnalysis::Transformations;

void getStr( size_t i, string &str ) {
	for ( size_t j = i ; j < i*i + 1; j++ ) {
		str.push_back( ( i*j ) % 256 );
	}
}

template < class Base = double > class Vector3 {
public:
	Base x, y, z;
};

int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	SummaryMultiplication<> a;
	FUNCTION_TYPE ft = a.GetType();
	UNREFERENCED_PARAMETER( ft );


	printf("Trololo Wololo...Testujeme...\n");

	//int barrier( 0 );
	//scanf_s( "%d", &barrier );

    return 0;
}

