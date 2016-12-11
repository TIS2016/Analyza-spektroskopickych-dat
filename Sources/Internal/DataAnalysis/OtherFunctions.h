#pragma once

enum SUMMARY_FUNCTION_TYPE {
	SFT_MULTIPLY = 0,
	SFT_DIVIDE = 1,
	SFT_UNKNOWN = 0xFFFFFFFF
};

class SummaryFunction : public IFunction {
public:
	SummaryFunction() : mType(SFT_UNKNOWN) {};
	SummaryFunction(__in SUMMARY_FUNCTION_TYPE type) : mType(type) {};

	virtual void Apply() {};

protected:
	SUMMARY_FUNCTION_TYPE mType;
};

class XTransform: public IFunction {
public:
	
	XTransform(__in double off, __in shared_ptr<IFunction> spScaleFunction ) : mOffset(off), mSpScaleFunction(spScaleFunction) {};

	virtual void Apply(__in double in, __out double &out) {
		mSpScaleFunction->Apply(in, out);
		out += mOffset;
	}

protected:
	double mOffset;
	shared_ptr<IFunction> mSpScaleFunction;
};

