#pragma once

#include "TransformationsLib.h"

namespace DataAnalysis { namespace Transformations {

	enum INPUT_TRANSFORMATION_STATE {
		ITS_NO_SAMPLES = 1,
		ITS_NO_TRANSFOMATIONS = 2,
		ITS_NO_MODEL = 4,
		ITS_TRANSFORM_DONE = 8,
		ITS_MODEL_DONE = 16
	};

	class InputTransformation {
	public:
		InputTransformation() {
			mState = ITS_NO_SAMPLES | ITS_NO_TRANSFOMATIONS | ITS_NO_MODEL;
			mXMax = numeric_limits<double>::lowest();
			mXMin = numeric_limits<double>::max();
		};

		void AddSample( __in const MeasurementSample &sample );
		void AddSample( __in const double x, __in const double y, __in const double uncertainty );

		void AddTransformation( __in const shared_ptr< IFunction<MeasurementSample> > spTransform );
		void AddTransformation( __in const TransformationHeader &transform );

		void CalculateTransformations( __in const size_t count, __out_ecount( count ) MeasurementSample *pOutput );
		void CalculateTransformations( __out Buffer<MeasurementSample> &output );
		void CalculateTransformations();

		void CalculateModel( __in const size_t count, __inout_ecount( count ) MeasurementSample *pOutput );
		void CalculateModel( __inout Buffer<MeasurementSample> &output );

		inline double GetXMin() const { return mXMin; };
		inline double GetXMax() const { return mXMax; };

	protected:

		int mState;

		vector< MeasurementSample > mInputSamples;
		vector< MeasurementSample > mOutputSamples;

		double mXMax;
		double mXMin;

		vector< shared_ptr< IFunction<MeasurementSample> > > mTranformations;
		vector< shared_ptr< IFunction<MeasurementSample> > > mModelTransformation;

	protected:

		inline bool FlagIsSet( __in const int flag ) const;

		inline void SetFlag( __in const int flag );

		inline void UnsetFlag( __in const int flag );

		inline void SetBoundValues();

	};
} }

