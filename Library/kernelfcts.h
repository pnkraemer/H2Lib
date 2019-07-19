/*
* A few kernel functions
*
* Author: Nicholas Krämer
*/


#ifndef KERNELFCTS_H
#define KERNELFCTS_H


#include "settings.h"


HEADER_PREFIX field
tps_kernel_1d(const real *xx, const real *yy, void *data);


HEADER_PREFIX field
exp_kernel_1d(const real *xx, const real *yy, void *data);



HEADER_PREFIX field
tps_kernel_2d(const real *xx, const real *yy, void *data);

HEADER_PREFIX field
exp_kernel_2d(const real *xx, const real *yy, void *data);


HEADER_PREFIX field
tps_kernel_s2(const real *xx, const real *yy, void *data);




#endif
