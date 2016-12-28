// DataAnalysisDebug.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "../Transformations/TransformationsLib.h"
using namespace DataAnalysis::Transformations;

// #include "../Utils/MathLib.h"
#include <iostream>

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

#include "../Transformations/HermitePolynomialProvider.h"

int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	uint cnstCnt = 3;
	int coefs[3] = { 1, 2, 3 };
	double coefs2[3] = { 1, 2, 3 };
	float coefs3[3] = { 1,2,3 };

	IFunction<double> *f1 = new LegendrePolynomialTransform<double>();
	IFunction<float> *f2 = new HermitePolynomialTransform<float>();
	IFunction<int> *f3 = new BasicPolynomialTransform<int>();

	f1->Initialize<LegendrePolynomialTransform<double>>( cnstCnt, coefs2 );
	f2->Initialize<HermitePolynomialTransform<float>>( cnstCnt, coefs3 );
	f3->Initialize<BasicPolynomialTransform<int>>( cnstCnt, coefs );

	// double f1o = 0;
	float f2o = 0;
	int f3o = 0;

	int in = 5;
	// f1->Apply( in, f1o );
	f2->Apply( (float)in, f2o );
	f3->Apply( in, f3o );
	
	// f3->Initialize<BasicPolynomialTransform<int>>( 3, coefs );

	/*
	B x;
	double y = x.foo( 1 );
	y = x.goo(  2.0, 5 );
	UNREFERENCED_PARAMETER( y );
	*/

	/*
	vector<double> b { 0, 1, 2, 3, 4 };
	vector<double> c( 5, 0 );
	IFunction<double> *a = new BasicPolynomial<>();

	vector<double>d { 9, 1 };
	a->Initialize<BasicPolynomial<>>( d );
	a->ApplyOnData( b.size(), b.data(), c.data() );
	
	delete a;
	*/
	
	// shared_ptr<IFunction> GetFunction( .... ) {....}; <- returns uninitialized function

	// size_t cnt = 3;
	// double pts[6] = { 4280, 2.8, 4288, 2.6, 4295, 2.8 };

	/*
	 size_t cnt = 7;
	 double pts[14] = { -3, -1, -2, -1, -1, -1, 0, 0, 1, 1, 2, 1, 3, 1 };

	shared_ptr< IFunction<double> > spSpline( new HermiteCubicSpline<>() );
	spSpline->Initialize<HermiteCubicSpline<>>( cnt, pts );
	
	double delta = .25;
	double start = -3; // 4275;
	double end = 3; // 4300;
	vector< pair<double, double> > f;

	double tmp = 0;
	for ( double x = start; x <= end; x += delta ) {
		spSpline->Apply( x, tmp );
		f.push_back( pair<double, double>( x, tmp ) );
	}

	FILE *file;
	fopen_s( &file, "F:\\School\\3.rocnik\\Winter\\TIS\\Analyza-spektroskopickych-dat\\Builds\\Debug\\testPlot1.txt", "w" );

	for ( auto it = f.begin(); it != f.end(); it++ ) {
		fprintf( file, "%f, %f\n", it->first, it->second );
	}
	

	fclose( file );

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

