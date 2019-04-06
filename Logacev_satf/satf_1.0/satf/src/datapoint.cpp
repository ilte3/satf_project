
#include <Rcpp.h>
#include <Rmath.h>
#include "satf.h"
using namespace Rcpp;

extern bool log_undefined_values;


inline bool valid_probability(double probability) {
    return (probability >= 0 && probability <= 1 && !std::isnan(probability));
}

CDataPoint::CDataPoint(int _n_responses_yes, int _n_responses, double _time, int _index, 
                       CDataPoint* last_datapoint) 
{
  _dbg_method("constructor", 1);
  _dbg((0, "last: %x", last_datapoint));

   n_responses_yes = _n_responses_yes;
   n_responses = _n_responses;
   time = _time;
   index = _index;
   mLastDatapoint = last_datapoint;
   relative_criterion = nan("");
}

void CDataPoint::UpdateLogLikGradient(std::vector<double>& cur_par, CCoefs& coefs, Rcpp::DoubleVector& dm_row) {
    _dbg_method("UpdateLogLikGradient", 1);
  
  // TODO: Only recompute when necessary.
  UpdateLogLik(cur_par);
  ComputeLogLikGradient(coefs, dm_row, cur_par[parameter_corr_mrsat]);
}

void CDataPoint::UpdateLogLik(std::vector<double>& cur_par)
{
  _dbg_method("UpdateLogLik", 1);

  ComputeDprime(cur_par[parameter_satf_asymptote], cur_par[parameter_satf_invrate],
               cur_par[parameter_satf_intercept]);
  ComputeCriterion(cur_par[parameter_bias_max], cur_par[parameter_bias_invrate],
                   cur_par[parameter_bias_intercept], cur_par[parameter_bias_min]);
  ComputeLogLik(cur_par[parameter_corr_mrsat]);
  _dbg((1, "dprime(asymptote=%f, invrate=%f, intercept=%f) = %f\n", cur_par[parameter_satf_asymptote], cur_par[parameter_satf_invrate],
               cur_par[parameter_satf_intercept], dprime.value));
  _dbg((1, "criterion(max=%f, invrate=%f, intercept=%f, min=%f) = %f\n", cur_par[parameter_bias_max], cur_par[parameter_bias_invrate],
                   cur_par[parameter_bias_intercept], cur_par[parameter_bias_min], criterion.value));
  _dbg((1, "mLogLik = %f\n", mLogLik));
}

void CDataPoint::Reset() {
  _dbg_method("Reset", 1);
  dprime.reset();
  criterion.reset();
  relative_criterion = nan("");
  mLogLik = nan("");
  
  mLogLikGradient = std::vector<double>();
  mRelativeCriterionGradient = std::vector<double>();
}


void CDataPoint::ComputeLogLik(double corr_mrsat, bool tolerate_imprecisions) 
{
  _dbg_method("ComputeLogLik", 1);

  assert(n_responses == 1 || corr_mrsat == 0.0 );
  assert(corr_mrsat == 0.0 || mLastDatapoint != NULL);

  if(std::isnan(corr_mrsat) || corr_mrsat == 0.0 || mLastDatapoint == NULL) {
  
      double pYes = _pnorm( -relative_criterion );
      mLogLik = _dbinom(n_responses_yes, n_responses, pYes, true );
      
      _dbg((0, "relative_criterion=%.2f, p_yes = %.2f", relative_criterion, pYes));
      _dbg((0, "log-lik=%.2f", mLogLik));

  } else {
      _dbg((0, "corr_mrsat=%.2f, relative_criterion=%.2f (response=%d)", corr_mrsat, relative_criterion, n_responses_yes ));
      _dbg((0, "last relative_criterion=%.2f (response=%d)", mLastDatapoint->relative_criterion, n_responses_yes));
      double p = pnorm_conditional(corr_mrsat, relative_criterion, mLastDatapoint->relative_criterion, 
                                n_responses_yes == 1, mLastDatapoint->n_responses_yes==1);
      mLogLik = log(p);
      
      _dbg((0, "log-lik=%.2f", mLogLik));
  }
}

void CDataPoint::ComputeLogLikGradient(CCoefs& coefs, Rcpp::DoubleVector& dm_row, double corr_mrsat) 
{
  _dbg_method("ComputeLogLikGradient", 1);
  assert(n_responses == 1);
  assert(mLastDatapoint != NULL);
  
  ComputeRelativeCriterionGradient(coefs, dm_row);
  
  bool independent = std::isnan(corr_mrsat) || corr_mrsat == 0.0 || mLastDatapoint == NULL;
  mLogLikGradient = std::vector<double>(coefs.size(), 0.0);
  for(int i=0; i < coefs.size(); i++) {
    double derivative = 0.0;
    if(independent) {
      derivative = IndependentLikDerivative(coefs, i, dm_row[i]);

    } else {
      derivative = DependentLikDerivative(coefs, i, dm_row[i], corr_mrsat);
    }
    mLogLikGradient[i] = (1/exp(mLogLik))*derivative;
    
  }
}

double ddbinom(int k, int n, double p) {
  return _choose(n, k)*(pow(p, k-1)*pow(1-p, n-k)*k - pow(p, k)*pow(1-p, n-k-1)*n + pow(p,k)*pow(1-p, n-k-1)*k);
}

double CDataPoint::IndependentLikDerivative(CCoefs& coefs, int coef_index, double dm_value) 
{
  _dbg_method("IndependentLikDerivative", 1);
  if(dm_value == 0.0)
    return 0.0;
    
  if(coefs.ParameterType(coef_index) == parameter_corr_mrsat)
    return 0.0;

  double pyes = _pnorm( -relative_criterion);
  double pyes_deriv = -_dnorm( -relative_criterion );
  double p_nresp_yes_deriv = ddbinom(n_responses_yes, n_responses, pyes);
  double derivative = p_nresp_yes_deriv*pyes_deriv*mRelativeCriterionGradient[coef_index];
  return derivative;
}

double CDataPoint::DependentLikDerivative(CCoefs& coefs, int coef_index, double dm_value, double corr_mrsat) 
{
  _dbg_method("DependentLikDerivative", 1);
  assert(n_responses == 1);
  assert(mLastDatapoint != NULL);
  if(dm_value == 0.0)
    return 0.0;
    
  double derivative;
  if(coefs.ParameterType(coef_index) == parameter_corr_mrsat) {
      assert(dm_value != 0);
      derivative = pnorm_conditional_derivative_rho(corr_mrsat, relative_criterion, mLastDatapoint->relative_criterion,
                                                    n_responses_yes == 1, mLastDatapoint->n_responses_yes==1);
      derivative = derivative*dm_value*coefs.TransformFn(coef_index, CCoefs::FnConstrainDerivative);
      
  } else {
      derivative = pnorm_conditional_derivative_criteria(corr_mrsat, relative_criterion, mLastDatapoint->relative_criterion, 
                                                mRelativeCriterionGradient[coef_index], mLastDatapoint->mRelativeCriterionGradient[coef_index],
                                                n_responses_yes == 1, mLastDatapoint->n_responses_yes==1);
  }
  return derivative;
}

void CDataPoint::ComputeRelativeCriterionGradient(CCoefs& coefs, Rcpp::DoubleVector& dm_row) {
  _dbg_method("ComputeRelativeCriterionGradient", 1);
  mRelativeCriterionGradient = std::vector<double>(coefs.size(), 0.0);
  for(int i=0; i < coefs.size(); i++) {
    mRelativeCriterionGradient[i] = RelativeCriterionDerivative(coefs, coefs.ParameterType(i), coefs.TransformFn(i, CCoefs::FnConstrainDerivative), dm_row[i]);
  }
}

double CDataPoint::RelativeCriterionDerivative(CCoefs& coefs, Parameter type, double transform_deriv, double dm_value)
{
  _dbg_method("RelativeCriterionDerivative", 1);
  
  if(dm_value == 0.0)
    return 0.0;
    
  double deriv;
  switch(type) 
  {
    case parameter_satf_asymptote:
        deriv = -dprime.PartialDerivativeByAsymptote();
        break;
          
    case parameter_satf_invrate:
        deriv = -dprime.PartialDerivativeByInvrate();
        break;
          
    case parameter_satf_intercept:
        deriv = -dprime.PartialDerivativeByIntercept();
        break;

    case parameter_bias_max:
        deriv = criterion.PartialDerivativeByAsymptote();
        break;
    
    case parameter_bias_invrate:
        deriv = criterion.PartialDerivativeByInvrate();
        break;
        
    case parameter_bias_intercept:
        deriv = criterion.PartialDerivativeByIntercept();
        break;
        
    case parameter_bias_min:
        deriv = criterion.PartialDerivativeByMin();
        break;
        
    case parameter_corr_mrsat:
    case parameter_invalid:
        deriv = nan("");
        break;
  }
  return deriv*dm_value*transform_deriv;
}

/*
Rcpp::DoubleVector CDataPoint::ComputeLogLikGradient(CCoefs& coefs)
{
    Rcpp::DoubleVector gradients;

    if(!Update(coefs)) {
      for(int i=0; i < coefs.size(); i++)
        gradients.push_back(R_NaN);
      return gradients;
    }

    double t = dprime.time;
    int datapoint_index = index;

    double criterion_tau = exp(-(t-criterion.intercept)/criterion.invrate);
    double criterion_at = (criterion.asymptote - criterion.min)*criterion_tau - criterion.asymptote + dprime.value; 
    double criterion_dnorm_at = _dnorm(criterion_at);
    double criterion_pnorm_at = _pnorm(criterion_at);
    double criterion_eta = criterion_pnorm_at*(-1+criterion_pnorm_at);

    double dprime_tau = exp(-(t-dprime.intercept)/dprime.invrate);
    double dprime_at = -criterion.value + dprime.asymptote - dprime.asymptote*dprime_tau;
    double dprime_dnorm_at = _dnorm(dprime_at);
    double dprime_pnorm_at = _pnorm(dprime_at);
    double dprime_eta = dprime_pnorm_at * (-1 + dprime_pnorm_at);

    std::vector<int>& coef_types = mDM->coef_types();

    for(int coef_idx = 0; coef_idx < coefs.size(); coef_idx++)
    {
        double gradient = 0;

        double w = (*mDM)(datapoint_index, coef_idx);
        double criterion_zeta = w * criterion_dnorm_at * criterion_tau;
        double dprime_zeta = w * dprime_dnorm_at * dprime.asymptote * dprime_tau;
        
        if(w != 0)
        {
          switch( coef_types[coef_idx] )
          {
              case parameter_satf_asymptote:
                if(t < dprime.intercept)
                  continue;

                gradient = -(dprime_dnorm_at * w) / dprime_eta;
                gradient = gradient * (n_responses_yes*(1- dprime_tau) 
                                    + (-1 + dprime_tau)*n_responses * dprime_pnorm_at );
                break;

              case parameter_satf_invrate:
                if(t < dprime.intercept)
                  continue;

                gradient = -dprime_zeta / ( dprime_eta * pow(dprime.invrate,2) );
                gradient = gradient * ((-t+dprime.intercept) * n_responses_yes 
                                       + n_responses * dprime_pnorm_at * (t-dprime.intercept));
                break;

              case parameter_satf_intercept:
                if(t < dprime.intercept)
                  continue;

                gradient = -dprime_zeta / ( dprime_eta * dprime.invrate );
                gradient = gradient * (-n_responses_yes + n_responses * dprime_pnorm_at);
                break;

              case parameter_bias_max:
                if(t < criterion.intercept)
                  continue;

                gradient = criterion_dnorm_at * w / criterion_eta;
                gradient = gradient * (n_responses_yes - n_responses*criterion_pnorm_at  + 
                                      (n_responses*criterion_pnorm_at - n_responses_yes)*criterion_tau );
                break;
 
            case parameter_bias_invrate:
              if(t < criterion.intercept)
                continue;
              {
                double t_minus_intercept = (t - criterion.intercept);
                double intercept_minus_t = (criterion.intercept - t);

                gradient = criterion_zeta / (criterion_eta * pow(criterion.invrate,2));
                gradient = gradient*(n_responses_yes*(criterion.asymptote*intercept_minus_t + criterion.min*t_minus_intercept)
                            + n_responses*criterion_pnorm_at*(criterion.asymptote*t_minus_intercept + criterion.min*intercept_minus_t));
              }
                break;

            case parameter_bias_intercept:
              if(t < criterion.intercept)
                continue;
               gradient = criterion_zeta / ( criterion_eta * criterion.invrate );
               gradient = gradient*(n_responses_yes*(-criterion.asymptote+criterion.min) + n_responses*criterion_pnorm_at*(criterion.asymptote-criterion.min));
              break;

            case parameter_bias_min:
              if(t > criterion.intercept) {
                  gradient = -criterion_zeta / criterion_eta;
                  gradient = gradient*( -n_responses_yes + n_responses*criterion_pnorm_at );

              } else  {
                  double dnorm_at = _dnorm(-criterion.value + dprime.value);
                  double pnorm_at = _pnorm(-criterion.value + dprime.value);
                  gradient = -( dnorm_at*w*(-n_responses_yes+n_responses*pnorm_at) );
                  gradient = gradient / ( pnorm_at*(-1+pnorm_at) );
              }
              break;

            case parameter_corr_mrsat:
                // No derivative for this parameter yet.
                gradient = nan("");
                break;

          };
          double x_Dtransform = coefs.TransformFn( coef_idx, CCoefs::FnConstrainDerivative );
          gradient = gradient * x_Dtransform;
          gradients.push_back( gradient );
      }
    }
    return gradients;
}
*/
