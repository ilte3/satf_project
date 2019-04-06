#ifndef __MATH_AUX_H__
#define __MATH_AUX_H__

#include <Rcpp.h>
#include <assert.h>

double pnorm_conditional(double rho, double relative_criterion, double last_relative_criterion, 
                         bool response_above_criterion, bool last_response_above_criterion);
double pnorm2d(double x_upper, double y_upper, double rho, double second_order=true);
double pnorm2d_derivative_by_ab(double x_upper, double y_upper, double rho, double x_upper_deriv, double y_upper_deriv, bool second_order=true);
double pnorm2d_derivative_by_rho(double x_upper, double y_upper, double rho, bool second_order=true);

inline double logodds2p(double lodds) { return( exp(lodds)/(1+exp(lodds)) ); }
inline double _dnorm(double x, double mu=0.0, double sigma=1.0, bool lg=false) { return ::Rf_dnorm4(x, mu, sigma, lg?1:0); }
inline double _pnorm(double x, double mu=0.0, double sigma=1.0, bool lt=true, bool lg=false) { return ::Rf_pnorm5(x, mu, sigma, lt?1:0, lg?1:0); }
inline double _ldnorm(double x, double mu=0.0, double sigma=1.0, bool lg=true) { return ::Rf_dnorm4(x, mu, sigma, lg?1:0); }
inline double _lpnorm(double x, double mu=0.0, double sigma=1.0, bool lt=true, bool lg=true) { return ::Rf_pnorm5(x, mu, sigma, lt?1:0, lg?1:0); }

inline double _dgamma(double x, double shp, double scl, bool lg=false) { return ::Rf_dgamma(x, shp, scl, lg?1:0); }
inline double _dbinom(double x, double n, double p, bool lg=false)     { return ::Rf_dbinom(x, n, p, lg?1:0); }
inline double _choose(int n, int k)     { return ::Rf_choose(n, k); }



double pnorm_conditional_derivative_criteria(double rho, double relative_criterion, double last_relative_criterion, 
                                              double relative_criterion_deriv, double last_relative_criterion_deriv, 
                                              bool response_above_criterion, bool last_response_above_criterion);
double pnorm_conditional_derivative_rho(double rho, double relative_criterion, double last_relative_criterion,
                                        bool response_above_criterion, bool last_response_above_criterion);


// negatively accelerated exponential
inline double NAE(double t, double asymptote, double invrate, double intercept) {
  if(asymptote == 0)      return 0.0;
  else if(t <= intercept) return 0.0;
  else if(invrate <= 0)   return nan("");
  else return asymptote*(1-exp(-(1/invrate)*(t-intercept)));
}

// shifted negatively accelerated exponential
inline double SNAE(double time, double asymptote, double invrate, double intercept, double min) {
	return NAE(time, asymptote-min, invrate, intercept) + min;
}

class SNAEPoint {
public:
  inline SNAEPoint() {
    reset();
  }
  
  inline SNAEPoint(double t, double fasymptote, double finvrate, double fintercept, double fmin=0.0) {
    Update(t, fasymptote, finvrate, fintercept, fmin);
  }
  inline void reset() {
    Update(nan(""), nan(""), nan(""), nan(""), nan(""));
    value = nan("");
  }

  inline void Update(double t, double fasymptote, double finvrate, double fintercept, double fmin=0.0) {
    time = t;
    asymptote = fasymptote;
    invrate = finvrate;
    intercept = fintercept;
    min = fmin;
    if(asymptote == fmin) value = fmin;
    else value = SNAE(time, asymptote, invrate, intercept, min);
  }
  
  inline double PartialDerivativeByAsymptote() {
    if(time < intercept) 
      return 0;
    else
      return 1-exp(-(time-intercept)/invrate);
  }
  inline double PartialDerivativeByInvrate() {
    if(time < intercept) 
      return 0;
    else {
      return -(asymptote-min)*(time-intercept)*exp(-(time-intercept)/invrate)/pow(invrate,2); //-(asymptote-min)*(time-intercept)*exp(-(time-intercept)/invrate)/(invrate*invrate);
    }
  }
  inline double PartialDerivativeByIntercept() {
    if(time < intercept) 
      return 0;
    else
      return -(asymptote-min)*exp(-(time-intercept)/invrate)/invrate;
  }
  inline double PartialDerivativeByMin() {
  if(time < intercept) 
      return 1;
    else
      return exp(-(time-intercept)/invrate);
  }
  
public:
  double time;
  double asymptote;
  double min;
  double invrate;
  double intercept;
  double value;
};


#endif // __MATH_AUX_H__
