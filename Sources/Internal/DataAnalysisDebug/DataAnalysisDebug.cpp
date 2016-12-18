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

class A {
public:
	virtual int foo( int a ) = 0;
	int foo( int a, int b ) {
		return a + b;
	}

	virtual int moo( int a, int b ) = 0;
	
	template <class T1, class T2> T1 goo( T1 a, T2 b ) {
		return moo( a, b );
	}
};

class B : public A {
public:
	// this hides foo with 2 arguments from A
	virtual int foo( int a ) {
		return 9 + a - a;
	}
	
	virtual int moo( int a, int b ) {
		return a * b;
	}
};


int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	
	/*
	B x;
	double y = x.foo( 1 );
	y = x.goo(  2.0, 5 );
	UNREFERENCED_PARAMETER( y );
	*/

	vector<double> b { 0, 1, 2, 3, 4 };
	vector<double> c( 5, 0 );
	IFunction<double> *a = new BasicPolynomial<>();

	vector<double>d { 9, 1 };
	a->Initialize<BasicPolynomial<>>( d );
	a->ApplyOnData( b.size(), b.data(), c.data() );
	
	delete a;
	
	// shared_ptr<IFunction> GetFunction( .... ) {....}; <- returns uninitialized function

	/*
	SummaryMultiplication<> a;
	FUNCTION_TYPE ft = a.GetType();
	UNREFERENCED_PARAMETER( ft );
	*/

	printf("Trololo Wololo...Testujeme...\n");

	//int barrier( 0 );
	//scanf_s( "%d", &barrier );

    return 0;
}

