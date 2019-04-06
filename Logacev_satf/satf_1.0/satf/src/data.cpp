
#include <Rcpp.h>
#include "satf.h"
using namespace Rcpp;


/******************************
 *   class CDataContainer     *
 ******************************/

CDataContainer::CDataContainer(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
                  Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
                  Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames):
                        mCoefConstraints(constraints, dm_ncoef)
{
    _dbg_method("constructor", 1);

    mObjectiveFunctionCalls = 0;

    // save unique rows of the design matrix with lists of
    // corresponding datapoint indices
    for(int row_idx=0; row_idx < dm.nrow(); row_idx++) {
      DoubleVector dm_row = dm(row_idx, _);
      _dbg((0, "LENGTH: %d", dm_row.size() ));
      AddDMRow(dm_row, row_idx);
    }
      
    std::string name_time = as< std::string >(cnames["time"]);
    DoubleVector time = as< DoubleVector >(data[ name_time ]);
  
    if( dv.containsElementNamed("response") ) {
      
        std::string name_response = as< std::string >(dv["response"]);
        IntegerVector response_yes = as< IntegerVector >(data[ name_response ]);
        CDataPoint* last_datapoint = NULL;
        for(int i=0; i < data.nrows(); i++) {
          CDataPoint* dpoint = new CDataPoint(response_yes[i], 1, time[i], i, last_datapoint);
          mDatapoints.push_back( dpoint );
          last_datapoint = dpoint;
        }
  
    } else if( dv.containsElementNamed("dprime") ) {
        // not implemented
  
    } else if( dv.containsElementNamed("n.responses") ) {
        std::string name_yes = as< std::string >(dv["n.responses.yes"]);
        std::string name_all = as< std::string >(dv["n.responses"]);
        IntegerVector n_responses_yes = as< IntegerVector >(data[ name_yes ]);
        IntegerVector n_responses = as< IntegerVector >(data[ name_all ]);
        CDataPoint* last_datapoint = NULL;
        for(int i=0; i < data.nrows(); i++) {
          CDataPoint* dpoint = new CDataPoint( n_responses_yes[i], n_responses[i], time[i], i, last_datapoint);
          mDatapoints.push_back( dpoint );
          last_datapoint = dpoint;
        }
    }
      
    mEnabled = std::vector<bool>(mUniqueDMRows.size(), true);
    
    mForceUpdate = true;
}

CDataContainer::~CDataContainer()
{
  for(size_t i=0; i < mDatapoints.size(); i++) {
      delete mDatapoints[i];
      mDatapoints[i] = NULL;
  }
}


int CDataContainer::FindDMRow(CDesignMatrixRow& row)
{
  _dbg_method("FindDMRow", 1);

  for(size_t j=0; j < mUniqueDMRows.size(); j++) {
    if(mUniqueDMRows[j] == row)
      return j;
  }
  return -1;
}

void CDataContainer::AddDMRow(DoubleVector& dm_row, int datapoint_index)
{
  _dbg_method("AddDMRow", 1);

  CDesignMatrixRow row(dm_row, mCoefConstraints);
  int idx = FindDMRow(row);
  if(idx == -1) {
    mUniqueDMRows.push_back(row);
    idx = mUniqueDMRows.size()-1;
  }
  mUniqueDMRows[idx].AddDatapointIndex(datapoint_index);
}


Rcpp::DoubleVector CDataContainer::ObjectiveFunction(DoubleVector& raw_coefs, bool by_row, bool force_update)
{
  _dbg_method("ObjectiveFunction", 1);

  if(mObjectiveFunctionCalls == 0)
    mForceUpdate = true;
  mForceUpdate = mForceUpdate || force_update;

  // TODO: Ensure on-demand updating works.
  mForceUpdate = true;

  // initialize coefficients
  CCoefs coefs =  CCoefs(raw_coefs, CCoefs::FnConstrain, false, mCoefConstraints); //, &mLastCoefs);

  int n_llvector_len = by_row ? mDatapoints.size() : 1;
  DoubleVector LLVector(n_llvector_len);
  for(int i=0; i < n_llvector_len; i++) {
      LLVector[i] = 0.0;
  }

  for(size_t i=0; i < mUniqueDMRows.size(); i++)
  {
    if(mEnabled[i] ) {
        _dbg((0, "updating parameter for enabled dm row %d/%d", i, mUniqueDMRows.size() )); 
        mUniqueDMRows[i].ObjectiveFunction(coefs, mForceUpdate, mDatapoints, &LLVector);
    } else {
        // Computation of log-likelihood on the disabled part of the dataset is meant to ensure that optimization on the selected part does not yield
        // estimates which are invalid for the disabled part. The result is discarded unless it's NaN.
        _dbg((0, "updating parameter for disabled dm row %d/%d", i, mUniqueDMRows.size() )); 
        bool LL_defined = mUniqueDMRows[i].ObjectiveFunction(coefs, mForceUpdate, mDatapoints, NULL);
        if(!by_row && !LL_defined) {
          LLVector[0] = R_NaN;
          return LLVector;
        }
    }
  }
  mLastCoefs = coefs;
  mForceUpdate = false;
  mObjectiveFunctionCalls++;
  
  _dbg((0, "returning")); 
  return LLVector;
}

Rcpp::NumericMatrix CDataContainer::ObjectiveFunctionGradient(DoubleVector& raw_coefs, bool by_row, bool force_update)
{
  _dbg_method("ObjectiveFunctionGradient", 1);

  if(mObjectiveFunctionCalls == 0)
    mForceUpdate = true;
  mForceUpdate = mForceUpdate || force_update;

  // TODO: Ensure on-demand updating works.
  mForceUpdate = true;
  
  // initialize coefficients
  CCoefs coefs =  CCoefs(raw_coefs, CCoefs::FnConstrain, false, mCoefConstraints); //, &mLastCoefs);

  int nrows = by_row ? mDatapoints.size() : 1;
  int ncols = coefs.size();
  NumericMatrix LLGradient = NumericMatrix(nrows, ncols);
  for(int i=0; i < nrows; i++) {
    for(int j=0; j < ncols; j++) {
      LLGradient(i, j) = 0.0;
    }
  }
  
  for(size_t i=0; i < mUniqueDMRows.size(); i++)
  {
    if(mEnabled[i] ) {
        _dbg((0, "updating parameter for enabled dm row %d/%d", i, mUniqueDMRows.size() )); 
        mUniqueDMRows[i].ObjectiveFunctionGradient(coefs, force_update, mDatapoints, &LLGradient);
    }
  }
  _dbg((0, "all gradients obtained")); 
  // TODO: the next line causes some sort of error, not sure why
  // mLastCoefs = coefs;

/*  mLastCoefs = coefs;
  mForceUpdate = false;
  mObjectiveFunctionCalls++;
*/
  return LLGradient;
}

/*
Rcpp::DoubleVector CDataContainer::ObjectiveFunctionGradient(DoubleVector& raw_coefs, bool by_row, bool tolerate_imprecisions)
{
  _dbg_method("ObjectiveFunctionGradient", 1);

  DoubleVector gradient(raw_coefs.size(), 0.0);
  
/*
// initialize coefficients
  CCoefs coefs =  CCoefs(raw_coefs, CCoefs::FnConstrain, false, mCoefConstraints);

  if( !mDM.HasParameter( parameter_corr_mrsat ) )  {
    for(size_t i=0; i < mDatapoints.size(); i++) {
      if( mEnabled[i] ) {
          gradient += mDatapoints[i].ComputeLogLikGradient(coefs);
      }
    }

  } else {
      DoubleVector invalid_gradient(raw_coefs.size(), R_NaN);
      return invalid_gradient;
  }
  return gradient;
  *//*
      DoubleVector invalid_gradient(raw_coefs.size(), R_NaN);
      return invalid_gradient;
}
*/
/*
void CDataContainer::DetermineZeroRows(LogicalVector& zero_columns, std::vector<bool>& row_selected, bool all)
{
    for(size_t i=0; i < row_selected.size(); i++) {
      if(all) row_selected[i] = true;
      else    row_selected[i] = false;
    }
    
    CharacterVector column_names = zero_columns.names();
    for(int i=0; i < zero_columns.size(); i++) 
    {
      std::string column_name = as<std::string>(column_names[i]);
      bool zero_column = zero_columns[i];
      int col_idx = mCoefConstraints.FindCoefIndex(column_name);

      for(size_t row_idx=0; row_idx < mUniqueDMRows.size(); row_idx++) 
      {
        double cur_value = mUniqueDMRows[row_idx][col_idx];
        if(all) {
            if(zero_column) row_selected[row_idx] = row_selected[row_idx] & (cur_value == 0.0);
            else            row_selected[row_idx] = row_selected[row_idx] & (cur_value != 0.0);
        } else {
            if(zero_column) row_selected[row_idx] = row_selected[row_idx] | (cur_value == 0.0);
            else            row_selected[row_idx] = row_selected[row_idx] | (cur_value != 0.0);
        }
      }
    }
}
*/

bool CDataContainer::SelectCoefSubset(Rcpp::CharacterVector& coefnames)
{
  // create a boolean vector in which all deselected coefficients are marked with true
  std::vector<bool> deselected_coefs(mCoefConstraints.size(), true);
  for(int i=0; i < coefnames.size(); i++) 
  {
    std::string cur_name = as<std::string>(coefnames[i]);
    int index = mCoefConstraints.FindCoefIndex( cur_name );
    if(index != -1)
      deselected_coefs[index] = false;
  }
  
  // disable all design matrix rows which use deselected coefficients
  for(size_t i=0; i < mUniqueDMRows.size(); i++) {
    if(mUniqueDMRows[i].CheckIfAnyFreeCoefsUsed(deselected_coefs))
      mEnabled[i] = false;
    else
      mEnabled[i] = true;
  }
  
  // force update of dprime, criterion and log-likelihood
  mForceUpdate = true;
  
  return true;
}

void CDataContainer::ResetSubset()
{
  // enable all design matrix rows
  for(size_t i=0; i < mEnabled.size(); i++)
    mEnabled[i] = true;
    
  // force update of dprime, criterion and log-likelihood
  mForceUpdate = true;
}

std::vector<int> CDataContainer::ReturnSelection()
{
  std::vector<int> selection;  
  for(size_t i=0; i < mUniqueDMRows.size(); i++) {
    if(mEnabled[i] ) {
      CDesignMatrixRow& cur_dmrow = mUniqueDMRows[i];
      std::vector<int>& datapoint_idxs = cur_dmrow.DatapointIndices();
      selection.insert(selection.end(), datapoint_idxs.begin(), datapoint_idxs.end());
    }
  }
  return selection;
}

