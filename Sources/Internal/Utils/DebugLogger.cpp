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

void DebugLogger::LogMessage(__in_z const char *pMessage) {
	if (mLogStream.is_open()) {
		const time_t currentTime = time(nullptr);
		tm locTime;
		localtime_s(&locTime, &currentTime);
		mLogStream << put_time(&locTime, "%Y-%m-%d %H:%M:%S") << ": " << pMessage << endl;
	}
}

void DebugLogger::LogMessage(__in const string& message) {
	LogMessage(message.data());
}

END_NAMESPACE