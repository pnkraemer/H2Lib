/*
* A few kernel functions
*
* Author: Nicholas Krämer
*/
#include "kernelfcts.h"

#include "basic.h"

field 
exp_kernel_1d(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;

    norm2 = REAL_SQR(xx[0] - yy[0]); //+ REAL_SQR(xx[1] - yy[1]);

    return REAL_EXP(-norm2);
}

field 
tps_kernel_1d(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;

    norm2 = REAL_SQR(xx[0] - yy[0]);// + REAL_SQR(xx[1] - yy[1]);
    if(norm2>0.0){
        return norm2* REAL_LOG(REAL_SQRT(norm2));
    }else{
        return 0.0;
    }
}


