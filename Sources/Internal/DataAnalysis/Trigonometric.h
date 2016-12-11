#pragma once

enum TRIGONOMETRIC_TYPE {
	TT_SIN = 0,
	TT_COS = 1,
	TT_UNKNOWN = 0xFFFFFFFF
};

class TrigonometricFunction : public IFunction {
public:
	TrigonometricFunction() : mType(TT_UNKNOWN) {};
	TrigonometricFunction(__in TRIGONOMETRIC_TYPE type) : mType(type) {};

protected:
	TRIGONOMETRIC_TYPE mType;
};

class SinFunction : public TrigonometricFunction {
public:
	SinFunction(__in double a, __in double b, __in double c) : TrigonometricFunction(TT_SIN), mA(a), mB(b), mC(c) {};

	virtual void Apply() {};

protected:
	double mA, mB, mC;
};
class CosFunction : public TrigonometricFunction {
public:
	CosFunction(__in double a, __in double b, __in double c) : TrigonometricFunction(TT_COS), mA(a), mB(b), mC(c) {};

	virtual void Apply() {};

protected:
	double mA, mB, mC;
};


