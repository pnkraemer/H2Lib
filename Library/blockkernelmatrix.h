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

struct _blockkernelmatrix {
    pamatrix kmat;    /* in compliance with kernelmatrix.h: only symmetric matrices */
    ph2matrix h2kmat;
    pamatrix pb;    /* see above, hence only one block*/
    uint rows;
    /*
     * TODO:
     * Should a blockkernelmatrix object have a reference
     * to the kernelmatrix object it stems from?
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
build_from_kernelmatrix_blockkernelmatrix(pkernelmatrix km, bool use_h2, pblock broot, pclusterbasis cb);



/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */
HEADER_PREFIX void
addeval_blockkernelmatrix_avector(field alpha, pcblockkernelmatrix bkm, pcavector src, pavector trg);



#endif
