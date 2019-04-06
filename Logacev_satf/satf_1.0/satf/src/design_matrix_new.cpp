
#include <Rcpp.h>
#include "satf.h"
using namespace Rcpp;


CDesignMatrixRow::CDesignMatrixRow(DoubleVector& elements, 
                                   CCoefConstraints& coef_constraints):
                                   mCurrentParameters(parameter_invalid, nan(""))
{
  _dbg_method("constructor", 1);

  mElements = elements;
  mCoefConstraints = &coef_constraints;
}

bool CDesignMatrixRow::UpdateParameters(CCoefs& coefs, bool force_update)
{
  _dbg_method("UpdateParameters", 1);
  bool params_changed = false;
  for(int p = 0; p < parameter_invalid; p++) {
    if( force_update || coefs.HasParameterChanged((Parameter)p, *this) ) {
      mCurrentParameters[p] = coefs.ComputeParameter((Parameter)p, *this);
      params_changed = true;
    }
  }
  return params_changed;
}

bool CDesignMatrixRow::CheckIfAnyFreeCoefsUsed(std::vector<bool>& coefs) {
  _dbg_method("CheckIfAnyFreeCoefsUsed", 1);

  for(size_t coef_idx=0; coef_idx < coefs.size(); coef_idx++) {
    if(coefs[coef_idx] && mElements[coef_idx] != 0.0 && !mCoefConstraints->IsCoefFixed(coef_idx))
      return true;
  }
  return false;
}

double CDesignMatrixRow::operator==(CDesignMatrixRow& other) const
{ 
  _dbg_method("operator==", 1);
  
  if(mElements.size() != other.mElements.size())
    return false;
    
  for(int i=0; i < other.mElements.size(); i++) {
    if( mElements[i] != other.mElements[i] )
      return false;
  }
  return true;
}

inline void SaveLogLik(DoubleVector* LLVector, double cur_LL, int idx) {
    if(LLVector == NULL)
      return;
    if(LLVector->size() == 1)
      (*LLVector)[0] += cur_LL;
    else
      (*LLVector)[idx] = cur_LL;
}


bool CDesignMatrixRow::ObjectiveFunction(CCoefs& coefs, bool force_update,
                                         std::vector<CDataPoint*>& datapoints, DoubleVector *LLVector)
{
    _dbg_method("ObjectiveFunction", 1);
  
    double sum_LL = 0.0;
    bool params_changed = UpdateParameters(coefs, force_update);
    _dbg((0, "parameters %s updated", (params_changed?"":"not") ));
    for(size_t j=0; j < mDatapointIndices.size(); j++)
    {
      int idx = mDatapointIndices[j];
      CDataPoint* cur_datapoint = datapoints[idx];
      _dbg((0, "updating datapoint %d (0x%x)", idx, cur_datapoint));
      if(params_changed) {
        cur_datapoint->UpdateLogLik(CurrentParameters());
        _dbg((0, "parameters updated"));
      } else {
        _dbg((0, "parameters up-to-date"));
      }
      double cur_LL = cur_datapoint->LogLik();
      sum_LL += cur_LL;
      SaveLogLik(LLVector, cur_LL, mDatapointIndices[j]);
   }
  _dbg((0, "done"));
   return !std::isnan(sum_LL);
}


inline void SaveLogLikGradient(NumericMatrix* LLGradient, std::vector<double>& gradient, int row_idx) {
    if(LLGradient == NULL)
      return;
    if(LLGradient->nrow() == 1) 
      row_idx = 0;
    for(size_t i=0; i < gradient.size(); i++) {
      (*LLGradient)(row_idx, i) += gradient[i];
    }
}

void CDesignMatrixRow::ObjectiveFunctionGradient(CCoefs& coefs, bool force_update,std::vector<CDataPoint*>& datapoints, 
                                                 NumericMatrix *LLGradient)
{
    _dbg_method("ObjectiveFunctionGradient", 1);
    bool params_changed = UpdateParameters(coefs, force_update);
    _dbg((0, "parameters %s updated", (params_changed?"":"not") ));
    for(size_t j=0; j < mDatapointIndices.size(); j++)
    {
      CDataPoint* cur_datapoint = datapoints[ mDatapointIndices[j] ];
      _dbg((0, "updating datapoint %d (0x%x)", mDatapointIndices[j], cur_datapoint));
      
      if(params_changed)
        cur_datapoint->UpdateLogLikGradient(CurrentParameters(), coefs, mElements);
      _dbg((0, "gradient updated"));
      std::vector<double>& gradient = cur_datapoint->LogLikGradient();
      SaveLogLikGradient(LLGradient, gradient, mDatapointIndices[j]);
      _dbg((0, "gradient saved"));
   }
     _dbg((0, "returning")); 

}

