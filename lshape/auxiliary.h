
#ifndef AUXILIARY_H
#define AUXILIARY_H

void
print_pointset(real **pts, uint npts, uint dim)
{
    printf("\n");
    for(uint i = 0; i < npts; i++){
        printf("(");
        printf("%.1e", pts[i][0]);
        for(uint j = 1; j < dim; j++){
            printf(", %.1e", pts[i][j]);
        }
        printf(")\n");
    }
    printf("\n");
}


void 
construct_random_2d(pkernelmatrix km, bool print_yes){

    uint i;
    uint npts;

    npts = km->points;
    
    for(i = 0; i < npts; i++){
        km->x[i][0] = FIELD_RAND() * 0.5 + 1.0;
        km->x[i][1] = FIELD_RAND() * 0.5 + 1.0;
        if(print_yes == 1)
            printf("(%.3f, %.3f), ", km->x[i][0], km->x[i][1] );   
    }
}



void 
construct_lattice_2d(pkernelmatrix km, bool print_yes){

    uint i;
    uint npts;

    npts = km->points;
    real genvec0 = 1.0;
    real genvec1 = 433461.0;
    
    for(i = 0; i < npts; i++){
        km->x[i][0] = fmod(genvec0 * i / (npts), 1.0);
        km->x[i][1] = fmod(genvec1 * i / (npts), 1.0);
        if(print_yes == 1)
            printf("(%.3f, %.3f), ", km->x[i][0], km->x[i][1] );   
    }
}


/*

void 
construct_lshape_2d(pkernelmatrix km, bool print_yes, uint L){

    uint i;
    uint npts;
    ptri2d *l_mesh;

    l_mesh = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (L + 1));
    l_mesh[0] = new_lshape_tri2d();
    for (i = 0; i < L; i++) {	
      l_mesh[i + 1] = refine_tri2d(l_mesh[i], NULL);
    }
    check_tri2d(l_mesh[L]);	
    p1_2d = new_tri2dp1(l_mesh[L]);	

    npts = km->points;
    real genvec0 = 1.0;
    real genvec1 = 433461.0;
    
    for(i = 0; i < npts; i++){
        km->x[i][0] = fmod(genvec0 * i / (npts), 1.0);
        km->x[i][1] = fmod(genvec1 * i / (npts), 1.0);
        if(print_yes == 1)
            printf("(%.3f, %.3f), ", km->x[i][0], km->x[i][1] );   
    }
    del_tri2d(l_mesh);
}

*/

#endif