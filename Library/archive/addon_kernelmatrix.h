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

HEADER_PREFIX void 
assemble_pblock(pckernelmatrix km, pamatrix mtrx);

HEADER_PREFIX void 
assemble_pblock_trans(pckernelmatrix km, pamatrix mtrx);

HEADER_PREFIX void 
assemble_big_kernelmatrix(pkernelmatrix km, pamatrix mat);

#endif