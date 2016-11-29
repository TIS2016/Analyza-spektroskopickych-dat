#include "stdafx.h"
#include "DebugLogger.h"
#include <iomanip>
#include <ctime>

BEGIN_NAMESPACE

DebugLogger::DebugLogger(__in_z const char *pLogFilePath ) {
	Initialize(pLogFilePath);
	if (mLogStream.is_open()) {
		mLogStream << "---------------- START OF LOG SESSION ----------------" << endl;
	}
}

DebugLogger::~DebugLogger() {
	mLogStream << "----------------- END OF LOG SESSION -----------------" << endl << endl;
	mLogStream.flush();
	mLogStream.close();
}

void DebugLogger::Initialize(__in_z const char *pPath) {
	mLogStream.open(pPath, ios_base::out | ios_base::app);
}

void DebugLogger::Log(__in_z const char *pMessage) {
	if (mLogStream.is_open()) {
		const time_t currentTime = time(nullptr);
		tm locTime;
		localtime_s(&locTime, &currentTime);
		mLogStream << mIndentString << put_time(&locTime, "%Y-%m-%d %H:%M:%S") << ": " << pMessage << endl;
	}
}

void DebugLogger::Log(__in const string& message) {
	Log(message.data());
}


void DebugLogger::LogFormatted( __in const char *pFormat, ... ) {
	char buff[2048];
	ZeroMemory(buff, 2048);

	va_list args;
	va_start(args, pFormat);
	vsnprintf_s(buff, 2047, pFormat, args);
	mLogStream << mIndentString << buff << endl;
	va_end(args);
}

void DebugLogger::BeginSection() {
	mIndentString.push_back('\t');
}

void DebugLogger::EndSection() {
	mIndentString.pop_back();
}

END_NAMESPACE
