
#include "satf_math.h"

/* Approximation of the conditional bivariate normal cdf, based on Albers&Kallenberg (1994), 
   "A simple approximation to the bivariate normal distribution with large correlation coefficient".
    Other approximations are cited in Albers & Kallenberg (1994):
    * Cox & Wermuth (1991), and
    * http://www.dtic.mil/dtic/tr/fulltext/u2/a125033.pdf (when correlation is low)
    * http://www.nasa.gov/centers/dryden/pdf/87929main_H-1120.pdf 
    * Continuous Multivariate Distributions, Models and Applications (Kotz et al.)
*/

namespace AlbersKallenberg1994 {
  
  // derivative of the normal distribution density function
  inline double _dddnorm( double x) {
    return _dnorm(x)*(x*x - 1);
  }
  inline double _ddnorm(double x) {
    return -(_dnorm(x)*x);
  }
  
  inline double f_theta(double rho) {
    return (sqrt(1 - (rho*rho) )/rho);
  }
  inline double f_c(double rho, double a, double b) {
    return (rho*a - b)/ sqrt(1-(rho*rho));
  }  
  inline double f_eta(double c) {
    return _pnorm(-c);
  }
  inline double f_c1(double c, double eta) {
    return (_dnorm(c) / eta);
  }
  inline double f_zeta(double a, double theta, double c, double c1) {
    return (a + theta*(c1 - c));
  }

  inline double _pnorm2d(double a, double b, double rho, bool smaller, bool second_order)
  {
      double theta = f_theta(rho);
      double c = f_c(rho, a, b);
      double eta = f_eta(c);
      double c1 = _dnorm(c) / eta;
      double zeta = f_zeta(a, theta, c, c1);

      double part1;
      if(smaller) part1 = _pnorm(b);
      else        part1 = 1-_pnorm(a);
      
      double part2 = eta*_pnorm(a);
      double part3 = - eta*_pnorm( zeta );
      double res = part1 + part2 + part3;
      if(second_order) {
        double part4a = -0.5*(theta*theta)*eta;
        double part4b = ( 1 + (c - c1)*c1 );
        double part4c = _ddnorm( zeta ); // _ddnorm(a);//zeta;//_ddnorm( zeta );
        double part4 = part4a*part4b*part4c;
        res += part4;
      }
      return res;
  }

inline double _pnorm2d_derivative_by_ab(double a, double b, double rho, 
                                          double a_deriv, double b_deriv, 
                                          bool smaller, bool second_order=true)
  {
     double rho_square = rho*rho;
     double theta = f_theta(rho);
     double alpha = (rho*a-b)/sqrt(1-rho_square);
//     double beta = a/sqrt(1-rho_square)+(rho*a-b)*rho/pow(1-rho_square, 1.5);
     double nu = _dnorm(alpha)/_pnorm(-alpha);
     double nu2 = _ddnorm(alpha)/_pnorm(-alpha);
     double gamma = nu-alpha;
     double omega = (rho*a_deriv - b_deriv) / sqrt(1-rho_square);
     double nu_square = nu*nu;
     double theta_square = theta*theta;
      double c = f_c(rho, a, b);
      double eta = f_eta(c);
      double c1 = _dnorm(c) / eta;
      double zeta = f_zeta(a, theta, c, c1);

double part1_deriv;
    if(smaller) part1_deriv = _dnorm(b)*b_deriv;
    else        part1_deriv = -_dnorm(a)*a_deriv;
//double mu=0.0, double sigma=1.0, bool lg=false) 
    double part2_deriv = -exp(_ldnorm(-alpha)+_lpnorm(a))*omega + exp(_lpnorm(-alpha)+_ldnorm(a))*a_deriv;
    double part3_deriv = -(-_dnorm(-alpha)*_pnorm(a+gamma*theta)*omega+_pnorm(-alpha)*_dnorm(a+gamma*theta)*(a_deriv+theta*omega*(nu2+nu_square-1)));
    double res = part1_deriv + part2_deriv + part3_deriv;
    if(second_order) {
      double part4a = -0.5*(theta*theta)*eta;
      double part4b = ( 1 + (c - c1)*c1 );
      double part4c = _ddnorm( zeta );
      
      double part4a_deriv = .5*theta_square*_dnorm(-alpha)*omega;
      double part4b1_deriv = omega*nu + alpha*omega*nu2 + alpha*omega*nu_square;
      double part4b2_deriv = -2*nu*nu2*omega - 2*nu_square*nu*omega;
      double part4b_deriv = part4b1_deriv + part4b2_deriv;
      double part4c_deriv = _dddnorm(zeta)*(a_deriv+theta*omega*(nu2+nu_square-1));
      double part4_deriv = part4a_deriv*part4b*part4c + part4a*part4b_deriv*part4c + part4a*part4b*part4c_deriv;
      res += part4_deriv;
    }
    return res;
  }

  // TODO: Clearly, some divergence between the derivate obtained numerical by compareDerivatives() and this derivative is due to the fact that a lot of floating point
  // arithmetic is done here. Fix it.
  inline double _pnorm2d_derivative_by_rho(double a, double b, double rho,
                                             bool smaller, bool second_order=true)
  {
     double rho_square = rho*rho;
     double theta = f_theta(rho);
     double alpha = (rho*a-b)/sqrt(1-rho_square);
     double beta = a/sqrt(1-rho_square)+(rho*a-b)*rho/pow(1-rho_square, 1.5);
     double nu = _dnorm(alpha)/_pnorm(-alpha);
     double nu2 = _ddnorm(alpha)/_pnorm(-alpha);
     double gamma = nu-alpha;
     
    double part2 = _dnorm(-(rho*a-b)/sqrt(1-rho*rho))*(-a/sqrt(1-rho*rho)-(rho*a-b)*rho/pow(1-rho*rho, 1.5))*_pnorm(a);
    double part3 = -(-_dnorm(-alpha)*beta*_pnorm(a+gamma*theta)+_pnorm(-alpha)*_dnorm(a+gamma*theta)*(-gamma/sqrt(1-rho_square)-theta*gamma/rho+theta*beta*(_ddnorm(alpha)/_pnorm(-alpha)+pow(_dnorm(alpha)/_pnorm(-alpha),2)-1)));
    double res = part2 + part3;
    if(second_order) {
        double theta_square = theta*theta;
        double nu_square = nu*nu;
        double part4 = -(-1.0*_pnorm(-alpha)*_ddnorm(a+theta*gamma)*(1-gamma*nu)/rho-theta_square*_pnorm(-alpha)*_ddnorm(a+theta*gamma)*(1-gamma*nu)/rho-.5*theta_square*_dnorm(-alpha)*beta*_ddnorm(a+theta*gamma)*(1-gamma*nu)+.5*theta_square*_pnorm(-alpha)*(_dddnorm(a+theta*gamma)*(-gamma/sqrt(1-rho_square)-theta*gamma/rho+theta*beta*(nu2+nu_square-1))*(1-gamma*nu)+_ddnorm(a+theta*gamma)*(beta*(1-nu2-nu_square)*nu-gamma*beta*nu2-nu_square*gamma*beta)));
        res = res + part4;
    }
    return res;
} 


  inline bool check_constraints(double a, double b, double rho) {
      assert( 0 < rho && rho < 1 );
      if( (rho*a-b) >= abs(rho*b - a)) return true;
      else                             return false;
  }
  
  double pnorm2d(double x_upper, double y_upper, double rho, bool second_order) {
      if( check_constraints(x_upper, y_upper, rho) ) {
        return _pnorm2d(x_upper, y_upper, rho, true, second_order);
      }
      if( check_constraints(y_upper, x_upper, rho) ){
        return _pnorm2d(y_upper, x_upper, rho, true, second_order);
      }
      if( check_constraints(-x_upper, -y_upper, rho) ) {
        return _pnorm2d(-x_upper, -y_upper, rho, false, second_order);
      }
      if( check_constraints(-y_upper, -x_upper, rho) ) {
        return _pnorm2d(-y_upper, -x_upper, rho, false, second_order);
      }
      return nan("");
  }
  
    double pnorm2d_derivative_by_rho(double x_upper, double y_upper, double rho, bool second_order) {
      if( check_constraints(x_upper, y_upper, rho) ) {
        return _pnorm2d_derivative_by_rho(x_upper, y_upper, rho, true, second_order);
      }
      if( check_constraints(y_upper, x_upper, rho) ){
        return _pnorm2d_derivative_by_rho(y_upper, x_upper, rho, true, second_order);
      }
      if( check_constraints(-x_upper, -y_upper, rho) ) {
        return _pnorm2d_derivative_by_rho(-x_upper, -y_upper, rho, false, second_order);
      }
      if( check_constraints(-y_upper, -x_upper, rho) ) {
        return _pnorm2d_derivative_by_rho(-y_upper, -x_upper, rho, false, second_order);
      }
      return nan("");
    }
    
    inline double pnorm2d_derivative_by_ab(double x_upper, double y_upper, double rho, double x_upper_deriv, double y_upper_deriv, bool second_order) {
      if( check_constraints(x_upper, y_upper, rho) ) {
        return _pnorm2d_derivative_by_ab(x_upper, y_upper, rho, x_upper_deriv, y_upper_deriv, true, second_order);
      }
      if( check_constraints(y_upper, x_upper, rho) ){
        return _pnorm2d_derivative_by_ab(y_upper, x_upper, rho, y_upper_deriv, x_upper_deriv, true, second_order);
      }
      if( check_constraints(-x_upper, -y_upper, rho) ) {
        return _pnorm2d_derivative_by_ab(-x_upper, -y_upper, rho, -x_upper_deriv, -y_upper_deriv, false, second_order);
      }
      if( check_constraints(-y_upper, -x_upper, rho) ) {
        return _pnorm2d_derivative_by_ab(-y_upper, -x_upper, rho, -y_upper_deriv, -x_upper_deriv, false, second_order);
      }
      return nan("");
    }
}

double pnorm2d(double x_upper, double y_upper, double rho, double second_order) {
    return AlbersKallenberg1994::pnorm2d(x_upper, y_upper, rho, second_order);
}

double pnorm2d_derivative_by_ab(double x_upper, double y_upper, double rho, double x_upper_deriv, double y_upper_deriv, bool second_order) {
  return AlbersKallenberg1994::pnorm2d_derivative_by_ab(x_upper, y_upper, rho, x_upper_deriv, y_upper_deriv, second_order);
}
  
double pnorm2d_derivative_by_rho(double x_upper, double y_upper, double rho, bool second_order) { 
  return AlbersKallenberg1994::pnorm2d_derivative_by_rho(x_upper, y_upper, rho, second_order);
}

 
double pnorm_conditional(double rho, double relative_criterion, double last_relative_criterion, 
                         bool response_above_criterion, bool last_response_above_criterion)
{
    if(last_response_above_criterion) {
      last_relative_criterion = -last_relative_criterion;
      relative_criterion = -relative_criterion;
    }
      
    double p_last = _pnorm(last_relative_criterion, 0.0, 1.0, true, false);
    double p_joint = pnorm2d(relative_criterion, last_relative_criterion, rho, true);
    
    if(response_above_criterion == last_response_above_criterion) {
      return p_joint/p_last;
    }
    return 1 - p_joint/p_last;
}

double pnorm_conditional_derivative_criteria(double rho, double relative_criterion, double last_relative_criterion, 
                                              double relative_criterion_deriv, double last_relative_criterion_deriv, 
                                              bool response_above_criterion, bool last_response_above_criterion)
{
    if(last_response_above_criterion) {
      last_relative_criterion = -last_relative_criterion;
      relative_criterion = -relative_criterion;
      relative_criterion_deriv = -relative_criterion_deriv;
      last_relative_criterion_deriv = -last_relative_criterion_deriv;
    }
    
    double p_last = _pnorm(last_relative_criterion, 0.0, 1.0, true, false);
    double p_joint = pnorm2d(relative_criterion, last_relative_criterion, rho, true);

    double deriv_last = _dnorm(last_relative_criterion)*last_relative_criterion_deriv;
    double deriv_joint = pnorm2d_derivative_by_ab(relative_criterion, last_relative_criterion, rho, relative_criterion_deriv, last_relative_criterion_deriv, true);
    
    if(response_above_criterion == last_response_above_criterion) {
      return (deriv_joint*p_last - deriv_last*p_joint)/(p_last*p_last);
    }
    // else
    return -deriv_joint/p_last+p_joint*deriv_last/(p_last*p_last);
}



double pnorm_conditional_derivative_rho(double rho, double relative_criterion, double last_relative_criterion, 
                                        bool response_above_criterion, bool last_response_above_criterion)
{
    
    if(last_response_above_criterion) {
      last_relative_criterion = -last_relative_criterion;
      relative_criterion = -relative_criterion;
    }
    
    double p_last = _pnorm(last_relative_criterion, 0.0, 1.0, true, false);
    double deriv_joint = pnorm2d_derivative_by_rho(relative_criterion, last_relative_criterion, rho, true);
    
    if(response_above_criterion == last_response_above_criterion) {
      return deriv_joint/p_last;
    }
    return - deriv_joint/p_last;
}


/*

P(X > x |Y > y) = P(X < -x |Y < -y)                     =     P(X < -x, Y < -y) / P(Y < -y)
P(X < x |Y > y) = P(X > -x |Y < -y) = P(X > -x |Y < -y) = 1 - P(X < -x, Y < -y) / P(Y < -y)
P(X > x |Y < y) = 1 - P(X < x |Y < y) = 1 - P(X < x, Y < y) / P(Y < y)
P(X < x |Y < y)                       =     P(X < x, Y < y) / P(Y < y)



*/

/*
  if( !valid_probability(pNo) ) {
    if(log_undefined_values) {
      printf("SATF WARNING: invalid conditional probability pNo=<%f>\n", pNo);
      printf("criterion=%f, dprime=%f, last_criterion=%f, last_dprime=%f\n", 
              criterion.value, dprime.value, last_datapoint.criterion.value, 
              last_datapoint.dprime.value);
      printf("corr.mrsat <%f>, relative_criterion <%.3f>, last_relative_criterion <%.3f>\n", corr_mrsat, relative_criterion, last_datapoint.relative_criterion);
      printf("last_response <%d>\n", last_datapoint.n_responses_yes);
    }
    return nan("");
  }
     

*/
