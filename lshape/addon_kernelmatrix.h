/*
* A few add-ons for the kernelmatrix.h module,
* mostly regarding conditionally positive definite kernels
*
* Author: Nicholas Krämer
*/

#ifndef ADDON_KERNELMATRIX
#define ADDON_KERNELMATRIX

#include "settings.h"
#include "kernelmatrix.h"

HEADER_PREFIX pamatrix 
assemble_pblock(pckernelmatrix km);

HEADER_PREFIX pamatrix 
assemble_pblock_trans(pckernelmatrix km);

HEADER_PREFIX void 
addeval_cond_kernelh2matrix(field alpha, ph2matrix h2km, pamatrix pb, pavector src, pavector trg);

HEADER_PREFIX void 
assemble_big_kernelmatrix(pkernelmatrix km, pamatrix mat);

#endif