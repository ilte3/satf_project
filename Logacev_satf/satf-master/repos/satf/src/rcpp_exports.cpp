
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "satf.h"
#include "satf_math.h"

#include "debug.h"

int global_dbg_level = 10;

static CDataContainer *zzz=NULL;

// [[Rcpp::export]]
LogicalVector rcpp_initialize_logLikFn(CharacterVector& dv, NumericMatrix& dm, IntegerVector& dm_ncoef, 
                                       NumericMatrix& constraints, DataFrame& data, CharacterVector& cnames)
{	
  _dbg_function("rcpp_initialize_logLikFn", 1);
    zzz = new CDataContainer(dv, dm, dm_ncoef, constraints, data, cnames);
    _dbg((0, "returning"));
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_update_constraints_logLikFn(NumericMatrix& constraints)
{  
  _dbg_function("rcpp_initialize_logLikFn", 1);
  if(zzz != NULL) {
    zzz->UpdateConstraints(constraints);
    return Rcpp::wrap(true);
  }
  return Rcpp::wrap(false);
  _dbg((0, "returning"));
}

// [[Rcpp::export]]
LogicalVector rcpp_deinitialize_logLikFn() {
  _dbg_function("rcpp_deinitialize_logLikFn", 1);
  _dbg((0, "pointer <0x%x>", zzz));
	delete zzz;
	zzz = NULL;
  _dbg((0, "returning"));
  return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_is_initialized_logLikFn() {
  _dbg_function("rcpp_is_initialized_logLikFn", 1);
  return Rcpp::wrap(zzz != NULL);
}

// [[Rcpp::export]]
CharacterVector rcpp_get_coef_names() {
  _dbg_function("rcpp_is_initialized_logLikFn", 1);
  if(zzz != NULL)
    return Rcpp::wrap(zzz->CoefNames());
  else
    return Rcpp::wrap(CharacterVector());
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, bool by_row=false, bool tolerate_imprecision=true, bool force_update=false) {
  _dbg_function("rcpp_compute_logLikFn", 1);
  if(zzz != NULL) {
    DoubleVector res = zzz->ObjectiveFunction(coefs, by_row, force_update);
    _dbg((0, "returning"));
    return Rcpp::wrap( res );
  } else {
    return Rcpp::wrap(DoubleVector());
  }
}


// [[Rcpp::export]]
NumericMatrix rcpp_compute_logLikFn_gradient(DoubleVector& coefs, bool by_row=false, bool tolerate_imprecision=true) {
  _dbg_function("rcpp_compute_logLikFn_gradient", 1);
  if(zzz != NULL) {
    NumericMatrix gradient = zzz->ObjectiveFunctionGradient(coefs, by_row, tolerate_imprecision);
    _dbg((0, "returning"));
    return Rcpp::wrap(gradient);
  }
  return Rcpp::wrap(DoubleVector());
}

// [[Rcpp::export]]
DoubleVector rcpp_constrain_coefs(DoubleVector& coefs) {
  _dbg_function("rcpp_constrain_coefs", 1);
  if(zzz != NULL) {
      DoubleVector res = zzz->ConstrainCoefs(coefs, true);
      _dbg((0, "returning"));
      return Rcpp::wrap(res);
  }
  return Rcpp::wrap(DoubleVector());
}

// [[Rcpp::export]]
DoubleVector rcpp_unconstrain_coefs(DoubleVector& coefs) {
  _dbg_function("rcpp_constrain_coefs", 1);
      if(zzz != NULL) {
      DoubleVector res = zzz->UnconstrainCoefs(coefs);
      _dbg((0, "returning"));
      return Rcpp::wrap(res);
  }
  return Rcpp::wrap(DoubleVector());
}

// [[Rcpp::export]]
bool rcpp_select_coef_subset(Rcpp::CharacterVector& coefnames) {
  _dbg_function("rcpp_select_coef_subset", 1);
  if(zzz != NULL) {
    bool res = zzz->SelectCoefSubset(coefnames);
    _dbg((0, "returning"));
    return Rcpp::wrap( res );
  }
  return Rcpp::wrap( false );
}

// [[Rcpp::export]]
void rcpp_reset_selection( ) {
  _dbg_function("rcpp_reset_selection", 1);
  if(zzz != NULL) {
    zzz->ResetSubset();
    _dbg((0, "returning"));
  }
}

// [[Rcpp::export]]
IntegerVector rcpp_return_selection( ) {
  _dbg_function("rcpp_reset_selection", 1);
  if(zzz != NULL) {
      std::vector<int> selection = zzz->ReturnSelection();
      // map to indices starting at 1, as used in R
      for(size_t i=0; i < selection.size(); i++)
        selection[i] = selection[i] + 1;
      _dbg((0, "returning"));
      return Rcpp::wrap( selection );
  }
  return Rcpp::wrap( IntegerVector() );
}

// [[Rcpp::export]]
void rcpp_set_coef_values(DoubleVector& values) {
  _dbg_function("rcpp_set_coef_values", 1);
  if(zzz != NULL) {
      zzz->SetCoefValues(values);
      _dbg((0, "returning"));
  }
}

// [[Rcpp::export]]
void rcpp_reset_coef_ranges( CharacterVector& names ) {
  _dbg_function("rcpp_reset_coef_values", 1);
  if(zzz != NULL) {
     zzz->ResetCoefRanges(names);
  _dbg((0, "returning"));
  }
}


// [[Rcpp::export]]
double rcpp_pnorm2d(double x_lower, double y_lower, double rho, bool second_order) {
    return pnorm2d(x_lower, y_lower, rho, second_order);
}

// [[Rcpp::export]]
double rcpp_pnorm2d_derivative_by_rho(double a, double b, double rho, bool second_order) {
    return pnorm2d_derivative_by_rho(a, b, rho, second_order);
}

// [[Rcpp::export]]
DoubleVector rcpp_correlate(IntegerVector& trial_id, DoubleVector& noise, DoubleVector& rho_vec) {
    double rho = rho_vec[0];
    for(int i=1; i < trial_id.length(); i++) {
      if(trial_id[i-1] == trial_id[i])
        noise[i] = noise[i-1]*rho + noise[i]*sqrt(1-pow(rho,2));
    }
  return Rcpp::wrap(noise);
}

/*
// [[Rcpp::export]]
NumericMatrix rcpp_get_dm() {
  return Rcpp::wrap(zzz->mDM.mDM);
}

// [[Rcpp::export]]
NumericVector rcpp_get_constraints_lower() {
  return Rcpp::wrap(zzz->mCoefConstraints.mCoefsLower);
}
// [[Rcpp::export]]
NumericVector rcpp_get_constraints_upper() {
  return Rcpp::wrap(zzz->mCoefConstraints.mCoefsUpper);
}
*/
