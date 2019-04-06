#include <Rcpp.h>

using namespace Rcpp;

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>
#include <vector>

#include "satf.h"
#include "satf_math.h"

bool log_undefined_values = true;

// TODO 1: Extend the transformation algorithm to include the possibility of stating that coefficient a
//         should not be smaller than another coefficient (b), or a should not be smaller than -1*b, etc. Without such
//         an option, we can end up with negative invrates, which is undesirable.

#define debug_level 10

_dbg_class_set(CCoefConstraints, "CCoefConstraints", 10);
_dbg_class_set(CCoefs, "CCoefs", 10);
_dbg_class_set(CDataPoint, "CDataPoint", 10);
_dbg_class_set(CDesignMatrixRow, "CDesignMatrixRow", 10);
_dbg_class_set(CDataContainer, "CDataContainer", 10);


int getMilliCount(){
  timeb tb;
	ftime(&tb);
	int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
	return nCount;
}

int getMilliSpan(int nTimeStart){
	int nSpan = getMilliCount() - nTimeStart;
	if(nSpan < 0)
		nSpan += 0x100000 * 1000;
	return nSpan;
}
/*
inline bool valid_dprime(double dprime, DoubleVector& LLVector, bool by_row) {
    if( isnan(dprime) ) {
        if(log_undefined_values)
          printf("SATF WARNING: dprime is undefined.\n");
        save_LL(LLVector, R_NaN, by_row);
        return false;
    }
    return true;
}
*/
