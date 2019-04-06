#ifndef __SATF_DATA_H__
#define __SATF_DATA_H__
#include <Rcpp.h>

//#define DEBUG
#include "debug.h"
#include "satf_math.h"


class CCoefs;

enum RVType {
  rv_binary = 1,
  rv_aggregate = 2,
  rv_dprime = 3
};

typedef enum EParameter {
  parameter_satf_asymptote = 0,
  parameter_satf_invrate  = 1,
  parameter_satf_intercept  = 2,
  parameter_bias_max  = 3,
  parameter_bias_invrate  = 4,
  parameter_bias_intercept  = 5,
  parameter_bias_min = 6,
  parameter_corr_mrsat  = 7,
  parameter_invalid = 8
} Parameter;



class CCoefConstraints;
class CDataPoint;

class CDesignMatrixRow
{
  public:
    CDesignMatrixRow(Rcpp::DoubleVector& elements, CCoefConstraints& coef_constraints);

    bool UpdateParameters(CCoefs& coefs, bool force_update);
    void AddDatapointIndex(int index) { mDatapointIndices.push_back(index); }

    std::vector<int>& DatapointIndices() { return mDatapointIndices; }
    std::vector<double>& CurrentParameters() { return mCurrentParameters; }
    
    bool ObjectiveFunction(CCoefs& coefs, bool force_update, std::vector<CDataPoint*>& datapoints, 
                            Rcpp::DoubleVector *LLVector);
    void ObjectiveFunctionGradient(CCoefs& coefs, bool force_update,std::vector<CDataPoint*>& datapoints, 
                                   Rcpp::NumericMatrix *LLGradient);

    bool CheckIfAnyFreeCoefsUsed(std::vector<bool>& coefs_set);

    double operator[](int i)  { return mElements[i]; }
    double operator==(CDesignMatrixRow& other) const;

  private:
    double ComputeParameter(CCoefs& coefs);
    
  private:
    Rcpp::DoubleVector mElements;
    std::vector<int> mDatapointIndices;
    std::vector<double> mCurrentParameters;
    CCoefConstraints *mCoefConstraints;

    _dbg_class_init;
};

class CDataPoint {
  public:
    CDataPoint(int _n_responses_yes, int _n_responses, double _time, 
               int _index, CDataPoint* last_datapoint);
  
    inline double LogLik() { return mLogLik; }
    inline std::vector<double>& LogLikGradient() { return mLogLikGradient; }

//  Rcpp::DoubleVector ComputeLogLikGradient(CCoefs& coefs);

    void UpdateLogLikGradient(std::vector<double>& cur_parameters, CCoefs& coefs, Rcpp::DoubleVector& dm_row);
    void UpdateLogLik(std::vector<double>& cur_parameters);
    void Reset();

private:
    void ComputeDprime(double satf_asymptote, double satf_invrate, double satf_intercept) {
      dprime.Update(time, satf_asymptote, satf_invrate, satf_intercept);
      relative_criterion = criterion.value - dprime.value;
    }
    // This is what Wickens (2001) calls lambda_center in equation (2.5) 
    void ComputeCriterion(double bias_max, double bias_invrate, double bias_intercept, double bias_min)  {
      criterion.Update(time, bias_max, bias_invrate, bias_intercept, bias_min);
      relative_criterion = criterion.value - dprime.value;
    }
    void ComputeLogLik(double corr_mrsat, bool tolerate_imprecisions=false);

    void ComputeLogLikGradient(CCoefs& coefs, Rcpp::DoubleVector& dm_row, double corr_mrsat);
    double IndependentLikDerivative(CCoefs& coefs, int coef_index, double dm_value);
    double DependentLikDerivative(CCoefs& coefs, int coef_index, double dm_value, double corr_mrsat);
    
    void ComputeRelativeCriterionGradient(CCoefs& coefs, Rcpp::DoubleVector& dm_row);
    double RelativeCriterionDerivative(CCoefs& coefs, Parameter type, double transform_deriv, double dm_value);
    
public:
    SNAEPoint dprime;
    SNAEPoint criterion;
    double relative_criterion;
    double mLogLik;
    
    int n_responses_yes;
    int n_responses;  
    double time;
    int index;
    
    std::vector<double> mRelativeCriterionGradient;
    std::vector<double> mLogLikGradient;
    CDataPoint* mLastDatapoint;

    _dbg_class_init;
};



class CCoefConstraints
{
  friend CCoefs;

  public:
    CCoefConstraints(Rcpp::NumericMatrix& coef_constraints, Rcpp::IntegerVector& dm_ncoef);
    ~CCoefConstraints();

    void UpdateConstraints(Rcpp::NumericMatrix& constraints);

    void SetCoefValues(Rcpp::DoubleVector& values);
    void ResetCoefRanges(Rcpp::CharacterVector& names);

    bool IsCoefFixed(int coef_index) { return (mCoefsLower[coef_index] == mCoefsUpper[coef_index]); }
    int  FindCoefIndex(std::string& name);

    int CoefIndex(std::string& column_name);
    bool HasParameter(Parameter parameter)  { return (mCoefIndexFirst[parameter] != -1); }
    int ParameterCoefIndexMin(Parameter parameter) { return mCoefIndexFirst[parameter]; }
    int ParameterCoefIndexMax(Parameter parameter) { return mCoefIndexLast[parameter]; }

    std::vector<int>& CoefTypes() { return mCoefTypes; }
    std::vector<std::string>& CoefNames() { return mCoefNames; }

    inline int size() { return mCoefNames.size(); } 

  private:
    bool SetCoefValue(std::string name, double value);
    bool ResetCoefRange(std::string name);

    Parameter StringToParameter(std::string& name);
    const char* ParameterToString(Parameter parameter);

  protected:
  //private:
    Rcpp::DoubleVector mCoefsLower;
    Rcpp::DoubleVector mCoefsUpper;
    std::vector<std::string> mCoefNames;

    Rcpp::DoubleVector mCoefsLowerOriginal;
    Rcpp::DoubleVector mCoefsUpperOriginal;
    
    std::vector<int> mCoefTypes;
    std::vector<int> mCoefIndexFirst;
    std::vector<int> mCoefIndexLast;

    _dbg_class_init;
};



// TODO: Rename 'Parameter' to 'ParameterType'
class CCoefs {
  public:
    typedef enum EFunction {
      FnConstrain,
      FnConstrainDerivative,
      FnUnconstrain
    } Function;

public:
    CCoefs();
    CCoefs(Rcpp::DoubleVector& coefs, Function fn, 
           bool use_names, CCoefConstraints& constraints,
           CCoefs* oldcoefs=NULL);

    Rcpp::DoubleVector Constrain(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector Unconstrain(Rcpp::DoubleVector& coefs);
    
    Rcpp::DoubleVector constrained() { return mConstrainedCoefs; }
    Rcpp::DoubleVector unconstrained() { return mUnconstrainedCoefs; }

    double ComputeParameter(Parameter parameter, CDesignMatrixRow& dm_row);
    bool HasParameterChanged(Parameter parameter, CDesignMatrixRow& dm_row);
    
    bool HasParameter(Parameter parameter) { return mConstraints->HasParameter(parameter); }

    double TransformFn(int coef_index, Function fn);
  
    Parameter ParameterType(int index) {return (Parameter)mConstraints->mCoefTypes[index];}
    Parameter ParameterType(std::string& name) {return mConstraints->StringToParameter(name);}

    inline int size() { return mConstrainedCoefs.size(); } 
    inline double operator[](int i) { return mConstrainedCoefs[i]; } 

private:
    double TransformFn(double x, double lower, double upper, Function fn);
    double TransformWithOneBoundary(double x, double lower, Function fn);
    double TransformWithTwoBoundaries(double x, double lower, double upper, Function fn);

private:
    Rcpp::DoubleVector mUnconstrainedCoefs;
    Rcpp::DoubleVector mConstrainedCoefs;
    std::vector<bool> mHasChanged;
    CCoefConstraints* mConstraints;

    _dbg_class_init;
};

class CDataContainer
{
  public:
    CDataContainer(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
              Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
              Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames);
    ~CDataContainer();

    void UpdateConstraints(Rcpp::NumericMatrix& constraints) { mCoefConstraints.UpdateConstraints(constraints); }

    Rcpp::DoubleVector ObjectiveFunction(Rcpp::DoubleVector& coefs, bool by_row=false, bool force_update=false);
    Rcpp::NumericMatrix ObjectiveFunctionGradient(Rcpp::DoubleVector& coefs, bool by_row=false, bool force_update=false);
    
    Rcpp::DoubleVector ConstrainCoefs(Rcpp::DoubleVector& raw_coefs, bool use_names)  {
      return CCoefs(raw_coefs, CCoefs::FnConstrain, use_names, mCoefConstraints).constrained();
    }
    Rcpp::DoubleVector UnconstrainCoefs(Rcpp::DoubleVector& raw_coefs)  {
      return CCoefs(raw_coefs, CCoefs::FnUnconstrain, false, mCoefConstraints).unconstrained();
    }

    bool SelectCoefSubset(Rcpp::CharacterVector& coefnames);
    void ResetSubset();
    std::vector<int> ReturnSelection();
    
    void SetCoefValues(Rcpp::DoubleVector& values)  { mCoefConstraints.SetCoefValues(values); }
    void ResetCoefRanges(Rcpp::CharacterVector& names)  { mCoefConstraints.ResetCoefRanges(names); }

    std::vector<std::string>& CoefNames() {return mCoefConstraints.CoefNames();}

  private:
    int  FindDMRow(CDesignMatrixRow& row);
    void AddDMRow(Rcpp::DoubleVector& row, int datapoint_index);

//  private:
  public:
    bool mForceUpdate;
  
    std::vector<CDesignMatrixRow> mUniqueDMRows;
    CCoefConstraints mCoefConstraints;
    std::vector<bool> mEnabled;

    CCoefs mLastCoefs;
    std::vector<CDataPoint*> mDatapoints;
    int mObjectiveFunctionCalls;

    _dbg_class_init;
};


#endif // __SATF_DATA_H__
