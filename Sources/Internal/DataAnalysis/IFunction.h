#pragma once
class IFunction {
public:
	
	virtual void Apply() = 0;
	virtual void Apply(__in double in, __out double &Out) = 0;
	virtual void Apply(__in size_t count, __in_ecount(count) double *pIn, __out_ecount(count) double *pOut) = 0;
	
};