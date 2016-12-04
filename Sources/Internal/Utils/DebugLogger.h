#pragma once
#include <fstream>

BEGIN_NAMESPACE

static const char* LOG_PATH_FILIP = "F:\\School\\3.rocnik\\Winter\\TIS\\Log01.txt";
static const char* MARTIN_LOG_PATH = "C:\\Users\\Martin\\Documents\\skola\\TIS_projekt\\Logy\\log.txt";

class DebugLogger
{
protected:
	ofstream mLogStream;
	string mIndentString;

protected:

	void Initialize(__in_z const char *pPath);

	void PrintRowHeader();

public:
	DebugLogger( __in_z const char *pLogFilePath );

	~DebugLogger();

	void Log(__in_z const char *pMessage);
	
	void Log(__in_z const string &message);

	void LogFormatted(__in_z const char *pFormat, ...);

	void BeginSection();

	void EndSection();

};

END_NAMESPACE