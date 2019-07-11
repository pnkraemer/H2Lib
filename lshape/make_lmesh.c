#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "matrixnorms.h"
#include "parameters.h"
#include "kernelmatrix.h"

#include "kernelfcts.h"
#include "addon_kernelmatrix.h"
#include "addon_krylovsolvers.h"

#include "../lshape/auxiliary.h"

#include "tri2d.h"      /* 2-dimensional mesh */












int
main(int argc, char **argv)
{
    init_h2lib(&argc, &argv);

    pstopwatch          sw;
    ptri2d              *gr_2d;
    uint                numpts, numref;
    uint                i, j;
    char                filename[250];
    FILE                *fptr;

    /* Parameters */
    sw = new_stopwatch();
    numref = atoi(argv[1]);     /* Number of refinements*/
    assert(numref<=8);           /* 8 refs are 200000 pts*/

    (void) printf("\nCreating L-shaped mesh,\n\t%u refinement(s)\n       ", numref);

    gr_2d = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (numref + 1));
    gr_2d[0] = new_lshape_tri2d();  
    for (i = 0; i < numref; i++)
        gr_2d[i + 1] = refine_tri2d(gr_2d[i], NULL);
    check_tri2d(gr_2d[numref]);     
    numpts =  gr_2d[numref]->vertices;
    printf("\t%u vertices\n", numpts);


    (void) printf("Saving mesh in file\n");
    sprintf(filename, "lshape/files/mesh/lmesh_N%d.txt", numpts);
    fptr = fopen(filename, "w"); 
    if (fptr == NULL) 
        printf("\n\nCOULD NOT OPEN FILE!\n\n"); 

    for (i=0; i<numpts; i++) 
        fprintf(fptr,"%lf %lf\n", gr_2d[numref]->x[i][0], gr_2d[numref]->x[i][1]); 
    fclose(fptr); 







    (void) printf("Cleaning up\n");
    for (i = 0; i <= numref; i++) {
      j = numref - i;
      del_tri2d(gr_2d[j]);
    }
    freemem(gr_2d);
    del_stopwatch(sw);

    uninit_h2lib();
    return 0;
}







