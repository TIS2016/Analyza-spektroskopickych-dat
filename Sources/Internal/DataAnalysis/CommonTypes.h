#pragma once

namespace DataAnalysis {
		 
		/*
			Container for spectroscopic measurement sample
				X -> X axis value
				Y -> Y axis value
				Deviation -> deviation (uncertainty) of measurement
				Model -> value of calculated model
		*/
		struct MeasurementSample {
			double X;
			double Y;
			double Deviation;
			double Model = 0;

			MeasurementSample() : X( 0 ), Y( 0 ), Deviation( 0 ), Model( 0 ) {};
			MeasurementSample( const double x, const double y, const double dev ) : X( x ), Y( y ), Deviation( dev ), Model( 0 ) {};
		};

}