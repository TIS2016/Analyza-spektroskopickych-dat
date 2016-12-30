#include "stdafx.h"
#include "InputParser.h"

namespace DataAnalysis { namespace Transformations {
	
  void Transformations::ConvertToInternal( __in const TransformationHeader &in, __out TransformationHeaderInternal &out ) {
  	out.transformationType = GetInternalFunctionType(in.name, 0);
  	out.functionParameters = in.functValues;
  
  	out.subFunctions.Allocate( in.subFunctions.Length() );
  	for (size_t i = 0; i < in.subFunctions.Length(); i++){
  		out.subFunctions[i] = GetInternalFunctionType( in.subFunctions[i], static_cast<int>(in.functValues[i][0]) );
  	}
  }

  shared_ptr<IFunction<double>> GetSubFunctionAndInitialize( __in const FUNCTION_TYPE type, __in const size_t paramCount, __in_ecount( paramCount ) const double *pParams, __in_z const char *pFunctName ) {
	  shared_ptr<IFunction<double>> spFunct = nullptr;

	  uint argCount = paramCount - 1;
	  const double *pArgs = pParams + 1;

	  switch ( type )
	  { 
	  case ( FT_POLY_BASIC ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new BasicPolynomialTransform<double>() );
		  _ASSERT( argCount > 0 );
		  if ( strncmp( "XScl", pFunctName, 4 ) == 0 ) {
			  spFunct->Initialize<BasicPolynomialTransform<double>>( paramCount, pParams ); // don't even ask, I'm about to cry
		  }
		  else {
			  spFunct->Initialize<BasicPolynomialTransform<double>>( argCount, pArgs );
		  }
		  break;
	  }
	  case ( FT_POLY_LEGENDRE ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new LegendrePolynomialTransform<double>() );
		  _ASSERT( argCount > 0 );
		  spFunct->Initialize<LegendrePolynomialTransform<double>>( argCount, pArgs );
		  break;
	  }
	  case ( FT_POLY_HERMITE ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new HermitePolynomialTransform<double>() );
		  _ASSERT( argCount > 0 );
		  spFunct->Initialize<HermitePolynomialTransform<double>>( argCount, pArgs );
		  break;
	  }
	  case ( FT_TRIG_SIN ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new SinFunction<double>() );
		  _ASSERT( argCount > 2 );
		  spFunct->Initialize<SinFunction<double>>( pArgs );
		  break;
	  }
	  case ( FT_TRIG_COS ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new SinFunction<double>() );
		  _ASSERT( argCount > 2 );
		  spFunct->Initialize<SinFunction<double>>( pArgs );
		  break;
	  }
	  case ( FT_SPLINE_CUBIC_HERMITE ):
	  {
		  spFunct = shared_ptr<IFunction<double>>( new HermiteCubicSpline<double>() );
		  _ASSERT( ( ( paramCount & 1) == 0 ) && ( paramCount > 1 ) );
		  uint cpCount = static_cast<uint>( paramCount >> 1 );
		  spFunct->Initialize<HermiteCubicSpline<double>>( cpCount, pParams );
		  break;
	  }
	  case ( FT_CONSTANT ):
		  spFunct = shared_ptr<IFunction<double>>( new ConstantFunction<double>() );
		  _ASSERT( paramCount == 1 );
		  spFunct->Initialize<ConstantFunction<double>>( *pParams );
		  break;
	  default:
		  break;
	  }

	  return spFunct;
  }

  shared_ptr<IFunction<double>> GetSubFunctionAndInitialize( __in const TransformationHeader &info, __in_z const char *pFunctName ) {
	  shared_ptr<IFunction<double>> spFunct = nullptr;

	  // search for given function in Header
	  for ( size_t i = 0; i < info.subFunctions.Length(); i++ ) {
		  string &name = info.subFunctions[i];
		  if ( name.compare( pFunctName ) == 0 ) {
			  // get internal representation of it (a.k.a. enum value)
			  FUNCTION_TYPE fType = GetInternalFunctionType( name, static_cast<int>( info.functValues[i][0] ) );
			  spFunct = GetSubFunctionAndInitialize( fType, info.functValues[i].Length(), info.functValues[i].Ptr(), pFunctName );
			  break;
		  }
	  }

	  return spFunct;
  }

  shared_ptr<IFunction<MeasurementSample>> GetSummaryFunction( __in const TransformationHeader &info ) {
	  shared_ptr<IFunction<MeasurementSample>> spFunct = nullptr;

	  // search for given function in Header
	  for ( size_t i = 0; i < info.subFunctions.Length(); i++ ) {
		  string &name = info.subFunctions[i];
		  if ( name.compare( "YTrg" ) == 0 ) {
			  // get internal representation of it (a.k.a. enum value)
			  FUNCTION_TYPE fType = GetInternalFunctionType( name, static_cast<int>( info.functValues[i][0] ) );

			  switch ( fType )
			  {
			  case FT_SUMOP_DIV:
				  spFunct = shared_ptr<IFunction<MeasurementSample>>( new SummaryDivision<double>() );
				  break;
			  case FT_SUMOP_MUL:
				  spFunct = shared_ptr<IFunction<MeasurementSample>>( new SummaryMultiplication<double>() );
				  break;
			  default:
				  break;
			  }

			  break;
		  }
	  }

	  return spFunct;
  }

  void CreateXTransform( __in const TransformationHeader &info, __inout const shared_ptr<IFunction<MeasurementSample>> spFunct ) {
	  /*
		1) get XOff sub-funct
		2) get XScl sub-funct
		3) Initialize XTransform with them
	  */
	  shared_ptr<IFunction<double>> spXOff = GetSubFunctionAndInitialize( info, "XOff" );
	  shared_ptr<IFunction<double>> spXScl = GetSubFunctionAndInitialize( info, "XScl" );

	  _ASSERT_EXPR( spXOff && spXScl, L"Nullptr in XTransform sub-functions" );

	  spFunct->Initialize<XTransform<>>( spXOff, spXScl );
  }

  void CreateYTransform( __in const TransformationHeader &info, __inout const shared_ptr<IFunction<MeasurementSample>> spFunct ) {
	  /*
		1) get YOff
		2) get YSpl, YTrg, YPol
		3) get YTyp (a.k.a. summary) 
		4) initialize YTyp with Yspl, YTrg, YPol
		5) initialize YTransform with YOff, YTyp
	  */

	  shared_ptr<IFunction<double>> spYOff = GetSubFunctionAndInitialize( info, "YOff" );
	  shared_ptr<IFunction<double>> spYSpl = GetSubFunctionAndInitialize( info, "YSpl" );
	  shared_ptr<IFunction<double>> spYTrg = GetSubFunctionAndInitialize( info, "YTrg" );
	  shared_ptr<IFunction<double>> spYPol = GetSubFunctionAndInitialize( info, "YPol" );

	  shared_ptr<IFunction<MeasurementSample>> spYSummary = GetSummaryFunction( info );

	  _ASSERT_EXPR( spYOff && spYSpl && spYTrg && spYPol && spYSummary, L"Nullptr in YTransform sub-functions" );
		   
	  spYSummary->Initialize<SummaryOperation<double>>( spYPol, spYTrg, spYSpl );
	  spFunct->Initialize<YTransform<double>>( spYOff, spYSummary );
  }


  shared_ptr<IFunction<MeasurementSample>> GetTransformation( __in const TransformationHeader &transformInfo ) {
	  shared_ptr<IFunction<MeasurementSample>> spFunct = nullptr;
	  
	  if ( transformInfo.name.compare( "XT" ) == 0 ) {
		  spFunct = shared_ptr<IFunction<MeasurementSample>>( new XTransform<>() );
		  CreateXTransform( transformInfo, spFunct );
	  }
	  else if ( transformInfo.name.compare( "YT" ) == 0 ) {
		  spFunct = shared_ptr<IFunction<MeasurementSample>>( new YTransform<>() );
		  CreateYTransform( transformInfo, spFunct );
	  }
	  else if ( transformInfo.name.compare( "BL" ) == 0 ) {
		  spFunct = shared_ptr<IFunction<MeasurementSample>>( new BaselineTransform<>() );
		  // CreateBaselineTransform( transformInfo, spFunct );
	  }
	  else if ( transformInfo.name.compare( "PK" ) == 0 ) {
		  // we'll see about this 0:)
	  }

	  return spFunct;
  }

} }