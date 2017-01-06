#include "stdafx.h"
#include "DebugLogger.h"
#include <iomanip>
#include <ctime>

namespace DataAnalysis { namespace Utils {

	DebugLogger::DebugLogger ( __in_z const char *pLogFilePath ) {
		Initialize ( pLogFilePath );
		if ( mLogStream.is_open () ) {
			mLogStream << "---------------- START OF LOG SESSION ----------------" << endl;
		}
	}

	DebugLogger::~DebugLogger () {
		mLogStream << "----------------- END OF LOG SESSION -----------------" << endl << endl;
		mLogStream.flush ();
		mLogStream.close ();
	}

	void DebugLogger::Initialize ( __in_z const char *pPath ) {
		mLogStream.open ( pPath, ios_base::out | ios_base::app );
	}

	void DebugLogger::PrintRowHeader () {
		const time_t currentTime = time ( nullptr );
		tm locTime;
		localtime_s ( &locTime, &currentTime );
		mLogStream << mIndentString << put_time ( &locTime, "%Y-%m-%d %H:%M:%S" ) << ": ";
	}

	void DebugLogger::Log ( __in_z const char *pMessage ) {
		if ( mLogStream.is_open () ) {
			PrintRowHeader ();
			mLogStream << pMessage << endl;
		}
	}

	void DebugLogger::Log ( __in const string& message ) {
		Log ( message.data () );
	}

	void DebugLogger::LogFormatted ( __in const char *pFormat, ... ) {
		if ( mLogStream.is_open () ) {
			char buff[2048];
			ZeroMemory ( buff, 2048 );

			va_list args;
			va_start ( args, pFormat );
			vsnprintf_s ( buff, 2047, pFormat, args );
			va_end ( args );

			PrintRowHeader ();
			mLogStream << buff << endl;
		}
	}

	void DebugLogger::BeginSection () {
		mIndentString.push_back ( '\t' );
	}

	void DebugLogger::EndSection () {
		if ( !mIndentString.empty () ) {
			mIndentString.pop_back ();
		}
	}

} }
