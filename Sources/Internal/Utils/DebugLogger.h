#pragma once
#include <fstream>

BEGIN_NAMESPACE

class DebugLogger
{
protected:
	ofstream mLogStream;

protected:

	void Initialize(__in_z const char *pPath);

public:
	DebugLogger( __in_z const char *pLogFilePath );

	~DebugLogger();

	void LogMessage(__in_z const char *pMessage);
	
	void LogMessage(__in_z const string &message);

};

END_NAMESPACE