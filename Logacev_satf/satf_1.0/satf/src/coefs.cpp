
#include <Rcpp.h>
#include "satf.h"
using namespace Rcpp;


/******************************
 *        class CCoefs        *
 ******************************/

CCoefs::CCoefs() 
{
}

CCoefs::CCoefs(DoubleVector& coefs, Function fn, bool use_names,
               CCoefConstraints& constraints, CCoefs* oldcoefs):
                                mHasChanged(coefs.size(), true)
{
  _dbg_method("constructor 1", 1);

  mConstraints = &constraints;

  if(fn == FnConstrain) {
    mUnconstrainedCoefs = coefs;
    mConstrainedCoefs = Constrain( coefs, use_names);

  } else {
    mConstrainedCoefs = coefs;
    mUnconstrainedCoefs = Unconstrain( coefs );
  }

  if(oldcoefs == NULL || oldcoefs->mUnconstrainedCoefs.size() != mUnconstrainedCoefs.size())
    return;

  for(int i = 0; i < mUnconstrainedCoefs.size(); i++)
  {
    if(oldcoefs->mUnconstrainedCoefs[i] == mUnconstrainedCoefs[i])
      mHasChanged[i] = false;
    else
      mHasChanged[i] = true;
  }
}


double CCoefs::TransformWithOneBoundary(double x, double lower, Function fn) {
  _dbg_method("TransformWithOneBoundary", 2);
  switch(fn) {
    case FnConstrain: // [-Inf, +Inf] -> [lower, +Inf]
      return pow(x, 2)+lower;

    case FnConstrainDerivative:
      return 2*x;

    case FnUnconstrain:  // [lower, +Inf] -> [-Inf, +Inf]
      return sqrt(x - lower);
  };
  return nan("");
}

double CCoefs::TransformWithTwoBoundaries(double x, double lower, double upper, Function fn) {
  _dbg_method("TransformWithTwoBoundaries", 2);
  
  switch(fn) {
    case FnConstrain: // [-Inf, +Inf] -> [0, 1] -> [lower, upper]
    {
        double y = exp(x)/(1+exp(x));
        double z = y*(upper-lower)+lower;
        return(z);
    }
    break;

    case FnConstrainDerivative:
    {
      return exp(x)*(upper-lower)/(1+exp(x)) - pow(exp(x),2)*(upper-lower)/pow(1+exp(x), 2);
    }
    break;
    
    case FnUnconstrain: // [lower, upper] -> [0, 1] -> [-Inf, +Inf]
    {
        if(x < lower || x > upper)
          return nan("");
      
        double y = (x-lower)/(upper-lower);
        double z = log(y/(1-y));
        
        if( std::isinf(z) ) {
          // transform of x is undefined because it is too close to the upper or lower boundary
          // so we will use 30 and -30 as the effective maxima
          z = ((x-lower) < (upper-x)) ? -30 : 30;
      }
      return(z);
    }
    break;
  };
  return nan("");
}


double CCoefs::TransformFn(double raw_coef, double lower, double upper, Function fn)
{
  _dbg_method("TransformFn 1", 1);

  if( lower == upper ) {
    switch(fn) {
      case FnConstrainDerivative:
        return 0;

      case FnConstrain:
      case FnUnconstrain:
        return lower;
    };
  
  } else if(lower == R_NegInf && upper == R_PosInf) {
    switch(fn) {
      case FnConstrainDerivative:
        return 1;

      case FnConstrain:
      case FnUnconstrain:
        return raw_coef;
    };

  } else if(upper == R_PosInf ) {
    return TransformWithOneBoundary(raw_coef, lower, fn);

  } else if(lower == R_NegInf) {
    return -1*TransformWithOneBoundary(-raw_coef, -upper, fn);
    
  } else {
    return TransformWithTwoBoundaries(raw_coef, lower, upper, fn);
    
  }
  return nan("");
}


double CCoefs::TransformFn(int coef_index, Function fn)
{
  _dbg_method("TransformFn 2", 1);
  switch(fn) {
    case FnConstrain:
    case FnConstrainDerivative:
        return TransformFn(mUnconstrainedCoefs[coef_index], mConstraints->mCoefsLower[coef_index], 
                         mConstraints->mCoefsUpper[coef_index], fn);

    case FnUnconstrain:
        return TransformFn(mConstrainedCoefs[coef_index], mConstraints->mCoefsLower[coef_index], 
                         mConstraints->mCoefsUpper[coef_index], fn);
  };
  return nan("");
}


Rcpp::DoubleVector CCoefs::Constrain(Rcpp::DoubleVector& raw_coefs, bool use_names) {
  _dbg_method("Constrain", 1);

  DoubleVector coefs; double val;

  for(int i=0; i < mConstraints->mCoefsLower.size(); i++)
  {
    val = TransformFn(raw_coefs[i], mConstraints->mCoefsLower[i], mConstraints->mCoefsUpper[i], FnConstrain);

    if(use_names) coefs.push_back( val, mConstraints->mCoefNames[i] );  
    else          coefs.push_back( val );  
  }
  return coefs;
}

DoubleVector CCoefs::Unconstrain(Rcpp::DoubleVector& raw_coefs)
{
  _dbg_method("Unconstrain", 1);

  DoubleVector coefs; double val;

  for(int i=0; i < raw_coefs.size(); i++) {
      val = TransformFn(raw_coefs[i], mConstraints->mCoefsLower[i], mConstraints->mCoefsUpper[i], FnUnconstrain);
      coefs.push_back( val );
  }
  return coefs;
}


// TODO: Update only when necessary.
double CCoefs::ComputeParameter(Parameter parameter, CDesignMatrixRow& dm_row)  
{
    _dbg_method("ComputeParameter", 1);

    if( !mConstraints->HasParameter(parameter) )
      return nan("");

    _dbg((0, "has parameter %s", mConstraints->ParameterToString(parameter) ));

    int min = mConstraints->ParameterCoefIndexMin(parameter);
    int max = mConstraints->ParameterCoefIndexMax(parameter);

    _dbg((0, "min %d, max %d\n", min, max));
    
    double param = 0.0;
    for(int i = min; i <= max; i++) 
         param += mConstrainedCoefs[i]*dm_row[i];

    _dbg((0, "returning %f", param));

    return param;
}


bool CCoefs::HasParameterChanged(Parameter parameter, CDesignMatrixRow& dm_row)
{
    _dbg_method("HasParameterChanged", 1);

/*
for(int i=0; i<mHasChanged.size(); i++)
  printf("%d, ", mHasChanged[i]?1:0);
printf("\n");

for(int i=0; i<mHasChanged.size(); i++)
  printf("%.2f, ", dm_row[i]);
printf("\n");
*/

    if( !mConstraints->HasParameter(parameter) ) {
      _dbg((0, "returning 0 because %s does not exist", mConstraints->ParameterToString(parameter)));
      return false;
    }

    int min = mConstraints->ParameterCoefIndexMin(parameter);
    int max = mConstraints->ParameterCoefIndexMax(parameter);
    
    bool changed = false;
    for(int i = min; i <= max; i++) 
    {
       if( std::isnan(dm_row[i]) || (dm_row[i] != 0.0 && mHasChanged[i]) ) {
          changed = true;
          break;
       }
    }
// printf("\n");
    _dbg((0, "returning %d for %s", changed, mConstraints->ParameterToString(parameter)));
    return changed;
}


/******************************
 *   class CCoefConstraints   *
 ******************************/

CCoefConstraints::CCoefConstraints(Rcpp::NumericMatrix& coef_constraints, Rcpp::IntegerVector& dm_ncoef):  
                                  mCoefTypes(coef_constraints.nrow(), -1),
                                  mCoefIndexFirst(parameter_invalid, -1),
                                  mCoefIndexLast(parameter_invalid, -1)
{
    _dbg_method("constructor", 1);
  
    UpdateConstraints(coef_constraints);
    // save coefficient names
    List dimnames(NumericMatrix(coef_constraints).attr("dimnames"));
    if (dimnames.size() > 0) {
        RObject names = dimnames[0];
        if (!names.isNULL()) {
            CharacterVector coef_names = CharacterVector(names);
            _dbg((0, "copying <%d> names", coef_names.size() ));
            for(int i=0; i < coef_names.size(); i++) {
                mCoefNames.push_back( as< std::string >(coef_names[i]) );
            }
        }
    }

    // create a mapping from parameter type to coefficient indices
    CharacterVector ncoef_names = dm_ncoef.names();
    int cur_coef_index = 0;
    for(int i=0; i < dm_ncoef.size(); i++) {
      std::string cur_name = as<std::string>(ncoef_names[i]); //mCoefNames[i];
      Parameter parameter = StringToParameter( cur_name );
      mCoefIndexFirst[parameter] = cur_coef_index;
      mCoefIndexLast[parameter] = mCoefIndexFirst[parameter] + dm_ncoef[i] - 1;
      _dbg((0, "parameter %s: indices [%d; %d]", cur_name.c_str(), mCoefIndexFirst[parameter], mCoefIndexLast[parameter]));
      cur_coef_index += dm_ncoef[i];
    }

    // create a mapping from coefficient index to parameter type
    for(size_t i=0; i < mCoefIndexFirst.size(); i++) {
        if(mCoefIndexFirst[i] < 0)
          continue;
          
        for(int j=mCoefIndexFirst[i]; j <= mCoefIndexLast[i]; j++) {
            mCoefTypes[j] = i;
        }
    }
}

void CCoefConstraints::UpdateConstraints(Rcpp::NumericMatrix& coef_constraints)
{
    _dbg_method("UpdateConstraints", 1);
  
    mCoefsLower = coef_constraints(_, 0);
    mCoefsUpper = coef_constraints(_, 1);

    // make a copy for SetCoefValues() / ResetCoefValues()
    mCoefsLowerOriginal = clone( mCoefsLower );
    mCoefsUpperOriginal = clone( mCoefsUpper );
}

CCoefConstraints::~CCoefConstraints() {
    _dbg_method("destructor", 1);
}

Parameter CCoefConstraints::StringToParameter(std::string& name)
{
  if(name == "asymptote") {
    return parameter_satf_asymptote;

  } else if(name == "invrate") {
    return parameter_satf_invrate;

  } else if(name == "intercept") {
    return parameter_satf_intercept;

  } else if(name == "bias.min") {
    return parameter_bias_min;

  } else if(name == "bias.max") {
    return parameter_bias_max;

  } else if(name == "bias.invrate") {
    return parameter_bias_invrate;
    
  } else if(name == "bias.intercept") {
    return parameter_bias_intercept;
    
  } else if(name == "corr.mrsat") {
    return parameter_corr_mrsat;
  }
  return parameter_invalid;
}

const char* CCoefConstraints::ParameterToString(Parameter parameter)
{

  switch(parameter) {
    case parameter_satf_asymptote:
      return "parameter_satf_asymptote";

    case parameter_satf_invrate:
      return "parameter_satf_invrate";

    case parameter_satf_intercept:
    return "parameter_satf_intercept";

    case parameter_bias_min:
      return "parameter_bias_min";

    case parameter_bias_max:
      return "parameter_bias_max";

    case parameter_bias_invrate:
      return "parameter_bias_invrate";
    
    case parameter_bias_intercept:
      return "parameter_bias_intercept";
    
    case parameter_corr_mrsat:
      return "parameter_corr_mrsat";

    case parameter_invalid:
      return "parameter_invalid";
  };
  return "unknown parameter";
}
int CCoefConstraints::CoefIndex(std::string& column_name) {
    for(size_t i=0; i < mCoefNames.size(); i++){
        if(mCoefNames[i] == column_name)
          return (int)i;
    }
    return -1;
}

bool CCoefConstraints::SetCoefValue(std::string name, double value) {
    int coef_idx = FindCoefIndex(name);
    if(coef_idx == -1)
      return false;

    mCoefsLower[coef_idx] = value;
    mCoefsUpper[coef_idx] = value;
    return true;
}

bool CCoefConstraints::ResetCoefRange(std::string name) {
    int coef_idx = FindCoefIndex(name);
    if(coef_idx == -1)
      return false;

    mCoefsLower[coef_idx] = mCoefsLowerOriginal[coef_idx];
    mCoefsUpper[coef_idx] = mCoefsUpperOriginal[coef_idx];
    return true;
}

int CCoefConstraints::FindCoefIndex(std::string& coef_name) {
  for(size_t i=0; i < mCoefNames.size(); i++) {
    if(mCoefNames[i] == coef_name)
      return i;
  }
  return -1;
}

void CCoefConstraints::SetCoefValues(DoubleVector& values) {
    CharacterVector names = values.names();
    for(int i=0; i < values.size(); i++) {
        std::string coef_name = as<std::string>(names[i]);
        SetCoefValue(coef_name, values[i]);
    }
}

void CCoefConstraints::ResetCoefRanges(CharacterVector& names) {
    for(int i=0; i < names.size(); i++) {
        std::string coef_name = as<std::string>(names[i]);
        ResetCoefRange(coef_name);
    }
}

