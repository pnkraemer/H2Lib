#ifndef BLOCKKERNELMATRIX_H
#define BLOCKKERNELMATRIX_H

typedef struct _blockkernelmatrix blockkernelmatrix;
typedef blockkernelmatrix *pblockkernelmatrix;
typedef const blockkernelmatrix *pcblockkernelmatrix;

#include "settings.h"
#include "h2matrix.h"
#include "kernelmatrix.h"
#include "clusterbasis.h"
#include "avector.h"
#include "h2arith.h"

struct _blockkernelmatrix {
    pamatrix kmat;      /* in compliance with kernelmatrix.h: only symmetric matrices */
    ph2matrix h2kmat;
    pamatrix pb;        /* see above, hence only one block */
    uint dof;           /* use "dof" because rows belongs to kmat */
    /*
     * TODO:
     * Should a blockkernelmatrix object have a reference
     * to the kernelmatrix object it stems from? (probably not)
     */
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

HEADER_PREFIX pblockkernelmatrix
new_blockkernelmatrix();

HEADER_PREFIX void
del_blockkernelmatrix(pblockkernelmatrix h2);

HEADER_PREFIX pblockkernelmatrix
build_from_kernelmatrix_full_blockkernelmatrix(pkernelmatrix km);

HEADER_PREFIX pblockkernelmatrix
build_from_kernelmatrix_h2_blockkernelmatrix(pkernelmatrix km, pblock broot, pclusterbasis cb);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

HEADER_PREFIX size_t
getsize_blockkernelmatrix(pcblockkernelmatrix bkm);


/* ------------------------------------------------------------
 * Compression
 * ------------------------------------------------------------ */

HEADER_PREFIX pblockkernelmatrix
compress_blockkernelmatrix(pcblockkernelmatrix bkm, real accuracy);


/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */
HEADER_PREFIX void
addeval_blockkernelmatrix_avector(field alpha, pcblockkernelmatrix bkm, pcavector src, pavector trg);

/* ------------------------------------------------------------
 * Conversion
 * ------------------------------------------------------------ */
HEADER_PREFIX void
convert_blockkernelmatrix_amatrix(pcblockkernelmatrix bkm, pamatrix mat);

/*
TODO: convert_blockkernelmatrix support for bkm if bkm is h2matrix
*/

#endif
