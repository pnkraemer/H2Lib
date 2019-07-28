/* ------------------------------------------------------------
 * TODO:
 * Does del_avector(top_subvec_trg) delete valuable coefficients?
 * ------------------------------------------------------------ */
#include "clusterbasis.h"
#include "blockkernelmatrix.h"
#include "basic.h"
#include "avector.h"
#include "h2compression.h"

/* ------------------------------------------------------------
 * Declarations of auxiliary functions
 * ------------------------------------------------------------ */
INLINE_PREFIX void
fill_pblock(pkernelmatrix km, pamatrix pb);

INLINE_PREFIX void
copy_blockkernel_topleft(pcblockkernelmatrix bkm, pamatrix mtrx);

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
build_from_kernelmatrix_full_blockkernelmatrix(pkernelmatrix km)  /* only one clusterbasis bc symmetric */
{
    uint points, dim;
    pblockkernelmatrix bkm;

    points = km->points;
    dim = km->dim;

    bkm = new_blockkernelmatrix();

    /* Kernel matrix */
    assert(bkm->h2kmat == NULL);
    bkm->kmat = new_amatrix(points, points);
    fillN_kernelmatrix(0, 0, km, bkm->kmat);

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

HEADER_PREFIX pblockkernelmatrix
build_from_kernelmatrix_h2_blockkernelmatrix(pkernelmatrix km, pblock broot, pclusterbasis cb)  /* only one clusterbasis bc symmetric */
{
    uint points, dim;
    pblockkernelmatrix bkm;

    points = km->points;
    dim = km->dim;

    bkm = new_blockkernelmatrix();

    /* H2 matrix */
    assert(bkm->kmat == NULL);
    assert(broot!=NULL);
    assert(cb!=NULL);
    bkm->h2kmat = build_from_block_h2matrix(broot, cb, cb);
    fill_h2matrix_kernelmatrix(km, bkm->h2kmat);

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
 * Statistics
 * ------------------------------------------------------------ */

HEADER_PREFIX size_t
getsize_blockkernelmatrix(pcblockkernelmatrix bkm)
{
    size_t sz;

    sz = (size_t) sizeof(blockkernelmatrix);

    if(bkm->h2kmat)
        sz += getsize_h2matrix(bkm->h2kmat);
    if(bkm->kmat)
        sz += getsize_amatrix(bkm->kmat);
    if(bkm->pb)
        sz += getsize_amatrix(bkm->pb);
    return sz;
}


/* ------------------------------------------------------------
 * Compression
 * ------------------------------------------------------------ */

HEADER_PREFIX pblockkernelmatrix
compress_blockkernelmatrix(pcblockkernelmatrix bkm, real accuracy)
{
    pblockkernelmatrix bkm_new;

    assert(bkm->h2kmat);
    assert(bkm->kmat == NULL);

    bkm_new = new_blockkernelmatrix();
    bkm_new->h2kmat = compress_h2matrix_h2matrix(bkm->h2kmat, false, false, 0, accuracy);

    if(bkm->pb)
    {
        bkm_new->pb = new_amatrix(bkm->pb->rows, bkm->pb->cols);
        copy_amatrix(0, bkm->pb, bkm_new->pb);
    }
    bkm_new->dof = bkm->dof;
    return bkm_new;
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

/* Auxiliary function which copies top left block of blockkernelmatrix object
into an amatrix object */
INLINE_PREFIX void
copy_blockkernel_topleft(pcblockkernelmatrix bkm, pamatrix mtrx)
{
    if(bkm->kmat){
        copy_amatrix(0, bkm->kmat, mtrx);
    }else if(bkm->h2kmat){
        const pamatrix fullmat = convert_h2matrix_amatrix(0, bkm->h2kmat);
        copy_amatrix(0, fullmat, mtrx);
        del_amatrix(fullmat);
    }
}

HEADER_PREFIX void
convert_blockkernelmatrix_amatrix(pcblockkernelmatrix bkm, pamatrix mat)
{
    uint points, pblen;
    pamatrix topleft, bottomleft, topright, bottomright;

    assert(bkm->dof == mat->rows);
    assert(mat->cols == mat->rows);

    /* If no polblock, then quick copy; otherwise block by block*/
    if(bkm->pb==NULL){
        copy_blockkernel_topleft(bkm, mat);
    }else{
        points = bkm->pb->rows;
        pblen = bkm->pb->cols;

        topleft = new_sub_amatrix(mat, points, 0, points, 0);
        bottomleft = new_sub_amatrix(mat, pblen, points, points, 0);
        topright = new_sub_amatrix(mat, points, 0, pblen, points);
        bottomright = new_sub_amatrix(mat, pblen, points, pblen, points);

        copy_blockkernel_topleft(bkm, topleft);
        copy_amatrix(1, bkm->pb, bottomleft);
        copy_amatrix(0, bkm->pb, topright);
        clear_amatrix(bottomright);

        del_amatrix(topleft);
        del_amatrix(bottomleft);
        del_amatrix(topright);
        del_amatrix(bottomright);
    }
}
