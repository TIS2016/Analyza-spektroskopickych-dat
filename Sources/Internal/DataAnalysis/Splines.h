#pragma once

enum SPLINE_TYPE {
	ST_CUBIC = 0,
	ST_UNKNOWN = 0xFFFFFFFF
};

// TODO: add abstract Spline class as a parent for all spline functions

struct CubicSplineParamPair {
	double a, b;
};

class CubicSpline : public IFunction {
public:
	CubicSpline(Buffer<CubicSplineParamPair> &splineParams) : mParameters(splineParams), mParameterCount(splineParams.Length()) {};

	virtual void Apply() {};

protected:
	size_t mParameterCount;
	Buffer<CubicSplineParamPair> mParameters;
};

