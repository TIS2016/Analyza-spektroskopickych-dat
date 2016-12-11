#pragma once

enum POLYNOMIAL_TYPE {
	PT_BASIC = 0,
	PT_LEGENDRE = 1,
	PT_CHEBYSHEV = 2,
	PT_UNKNOWN = 0xFFFFFFFF
};

class Polynomial : public IFunction {
public:
	Polynomial() : mType(PT_UNKNOWN) {};
	Polynomial(__in POLYNOMIAL_TYPE type) : mType(type) {};

protected:
	POLYNOMIAL_TYPE mType;
};

class BasicPolynomial : public Polynomial {
public:
	BasicPolynomial(Buffer<double> &vals) : Polynomial(PT_BASIC), mValues(vals), mDegree(vals.Length()) {};

	virtual void Apply() {};
protected:
	UINT mDegree;
	Buffer<double> mValues;
};

class LegendrePolynomial : public Polynomial {

};

class ChebyshevPolynomial : public Polynomial {

};