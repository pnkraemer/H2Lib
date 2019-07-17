/* ------------------------------------------------------------
 * TODO:
 * Does del_avector(top_subvec_trg) delete valuable coefficients?
 * ------------------------------------------------------------ */
#include "clusterbasis.h"
#include "blockkernelmatrix.h"
#include "basic.h"
#include "avector.h"


INLINE_PREFIX void
fill_pblock(pkernelmatrix km, pamatrix pb);

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

HEADER_PREFIX pblockkernelmatrix
new_blockkernelmatrix()
{
    pblockkernelmatrix bkm;
    bkm = allocmem(sizeof(blockkernelmatrix));
    bkm->kmat = NULL;
    bkm->h2kmat = NULL;
    bkm->pb = NULL;
    bkm->rows = 0;
    return bkm;
}

HEADER_PREFIX void
del_blockkernelmatrix(pblockkernelmatrix bkm)
{
    if (bkm->kmat) {
        del_amatrix(bkm->kmat);
    }

    if (bkm->h2kmat) {
        del_h2matrix(bkm->h2kmat);
    }

    if (bkm->pb) {
        del_amatrix(bkm->pb);
    }

    freemem(bkm);
}


INLINE_PREFIX void
fill_pblock(pkernelmatrix km, pamatrix pb)
{
    uint points, dim;
    uint rows, cols;
    uint i, j;

    points = km->points;
    dim = km->dim;
    rows = pb->rows;
    cols = pb->cols;
    assert(points == rows);
    assert(dim + 1 == cols);

    for(i = 0; i < points; i++){
        setentry_amatrix(pb, i, 0, 1.0);
        for(j = 0; j < dim; j++){
            setentry_amatrix(pb, i, 1 + j, km->x[i][j]);
        }
    }
}


HEADER_PREFIX pblockkernelmatrix
build_from_kernelmatrix_blockkernelmatrix(pkernelmatrix km, bool use_h2, pblock broot, pclusterbasis cb)  /* only one clusterbasis bc symmetric */
{
    uint points, dim;
    pblockkernelmatrix bkm;

    points = km->points;
    dim = km->dim;

    bkm = new_blockkernelmatrix();

    if(use_h2 == true){
        assert(cb!=NULL);
        bkm->h2kmat = build_from_block_h2matrix(broot, cb, cb);
        fill_h2matrix_kernelmatrix(km, bkm->h2kmat);
        printf("\n\nTODO: fill bkm->h2kmat with proper values!\n\n");
//        assert(false);
    }else{
        bkm->kmat = new_amatrix(points, points);
        fillN_kernelmatrix(0, 0, km, bkm->kmat);
    }
    if((km->cpos)>0){
        bkm->pb = new_amatrix(points, 1 + dim);
        fill_pblock(km, bkm->pb);
        bkm->rows = points + 1 + dim;
    }else{
        bkm->rows = points;
    }
    return bkm;
}



/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

HEADER_PREFIX void
addeval_blockkernelmatrix_avector(field alpha, pcblockkernelmatrix bkm, pcavector src, pavector trg)
{
    pavector    copy_src, copy_trg;
    pavector    top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint        pts, pbsize;


    if(bkm->pb){
        assert(src->dim == bkm->rows);

        copy_src = new_avector(src->dim);
        copy_avector(src, copy_src);
        copy_trg = new_avector(trg->dim);
        copy_avector(trg, copy_trg);

        pts = bkm->pb->rows;
        pbsize = bkm->pb->cols;

        top_subvec_src = new_sub_avector(copy_src, pts, 0);
        bottom_subvec_src = new_sub_avector(copy_src, pbsize, pts);
        top_subvec_trg = new_sub_avector(copy_trg, pts, 0);
        bottom_subvec_trg = new_sub_avector(copy_trg, pbsize, pts);

        if(bkm->h2kmat){
//            addeval_h2matrix_avector(alpha, bkm->h2kmat, top_subvec_src, top_subvec_trg);
        }
        if(bkm->kmat){
            addeval_amatrix_avector(alpha, bkm->kmat, top_subvec_src, top_subvec_trg);
        }
        addevaltrans_amatrix_avector(alpha, bkm->pb, top_subvec_src, bottom_subvec_trg);
        copy_avector(copy_trg, trg);

        del_avector(top_subvec_src);
        del_avector(bottom_subvec_src);
        del_avector(top_subvec_trg);
        del_avector(bottom_subvec_trg);
        del_avector(copy_src);
        del_avector(copy_trg);
    }
    else if(bkm->h2kmat){
        assert(src->dim == bkm->rows);
        addeval_h2matrix_avector(alpha, bkm->h2kmat, src, trg);
    }else if(bkm->kmat){
        assert(src->dim == bkm->rows);
        addeval_amatrix_avector(alpha, bkm->kmat, src, trg);
    }
}
