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
    bkm->dof = 0;
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

    /* H matrix or kernel matrix */
    if(use_h2 == true){
        assert(cb!=NULL);
        bkm->h2kmat = build_from_block_h2matrix(broot, cb, cb);
        fill_h2matrix_kernelmatrix(km, bkm->h2kmat);
    }else{
        bkm->kmat = new_amatrix(points, points);
        fillN_kernelmatrix(0, 0, km, bkm->kmat);
    }

    /* Polynomial part or not */
    if((km->cpos)>0){
        bkm->pb = new_amatrix(points, 1 + dim);
        fill_pblock(km, bkm->pb);
        bkm->dof = points + 1 + dim;
    }else{
        bkm->dof = points;
    }
    return bkm;
}



/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */
/*
HEADER_PREFIX void
addeval_blockkernelmatrix_avector(field alpha, pcblockkernelmatrix bkm, pcavector src, pavector trg)
{
    pavector    copy_src;
    pavector    top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint        pts, pbsize;


    if(bkm->pb){
        assert(src->dim == bkm->dof);

        // Copy source bc we need pavector and not pcavector
        copy_src = new_avector(src->dim);
        copy_avector(src, copy_src);

        pts = bkm->pb->rows;
        pbsize = bkm->pb->cols;

        top_subvec_src = new_sub_avector(copy_src, pts, 0);
        bottom_subvec_src = new_sub_avector(copy_src, pbsize, pts);
        top_subvec_trg = new_sub_avector(trg, pts, 0);
        bottom_subvec_trg = new_sub_avector(trg, pbsize, pts);

        if(bkm->h2kmat){
            addeval_h2matrix_avector(alpha, bkm->h2kmat, top_subvec_src, top_subvec_trg);
        }
        if(bkm->kmat){
            addeval_amatrix_avector(alpha, bkm->kmat, top_subvec_src, top_subvec_trg);
        }
        addevaltrans_amatrix_avector(alpha, bkm->pb, top_subvec_src, bottom_subvec_trg);


        del_avector(top_subvec_src);
        del_avector(bottom_subvec_src);
        del_avector(top_subvec_trg);
        del_avector(bottom_subvec_trg);
        del_avector(copy_src);
    }
    else if(bkm->h2kmat){
        assert(src->dim == bkm->dof);
        addeval_h2matrix_avector(alpha, bkm->h2kmat, src, trg);
    }else if(bkm->kmat){
        assert(src->dim == bkm->dof);
        addeval_amatrix_avector(alpha, bkm->kmat, src, trg);
    }
}
*/

HEADER_PREFIX void
addeval_blockkernelmatrix_avector(field alpha, pcblockkernelmatrix bkm, pcavector src, pavector trg)
{
    pavector    copy_src, copy_trg;
    pavector    top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint        pts, pbsize;


    if(bkm->pb){
        assert(src->dim == bkm->dof);

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
            addeval_h2matrix_avector(alpha, bkm->h2kmat, top_subvec_src, top_subvec_trg);
            addeval_amatrix_avector(alpha, bkm->pb, bottom_subvec_src, top_subvec_trg);
        }
        if(bkm->kmat){
            addeval_amatrix_avector(alpha, bkm->kmat, top_subvec_src, top_subvec_trg);
            addeval_amatrix_avector(alpha, bkm->pb, bottom_subvec_src, top_subvec_trg);
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
        assert(src->dim == bkm->dof);
        addeval_h2matrix_avector(alpha, bkm->h2kmat, src, trg);
    }else if(bkm->kmat){
        assert(src->dim == bkm->dof);
        addeval_amatrix_avector(alpha, bkm->kmat, src, trg);
    }
}


/* ------------------------------------------------------------
 * Conversion
 * ------------------------------------------------------------ */
HEADER_PREFIX void
convert_blockkernelmatrix_amatrix(pcblockkernelmatrix bkm, pamatrix mat)
{
    assert(bkm->h2kmat==NULL);  /*So far no h2matrix support (TODO!) */
    uint points, pblen;
    pamatrix topleft, bottomleft, topright, bottomright;

    assert(bkm->dof == mat->rows);
    assert(mat->cols == mat->rows);

    /* If no polblock, then quick copy; otherwise block by block*/
    if(bkm->pb==NULL){
        copy_amatrix(0, bkm->kmat, mat);
    }else{
        points = bkm->kmat->rows;
        pblen = bkm->pb->cols;
        topleft = new_sub_amatrix(mat, points, 0, points, 0);
        bottomleft = new_sub_amatrix(mat, pblen, points, points, 0);
        topright = new_sub_amatrix(mat, points, 0, pblen, points);
        bottomright = new_sub_amatrix(mat, pblen, points, pblen, points);

        copy_amatrix(0, bkm->kmat, topleft);
        copy_amatrix(1, bkm->pb, bottomleft);
        copy_amatrix(0, bkm->pb, topright);

        clear_amatrix(bottomright);
        del_amatrix(topleft);
        del_amatrix(bottomleft);
        del_amatrix(topright);
        del_amatrix(bottomright);

    }
}
