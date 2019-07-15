#include <stdio.h>
#include <string.h>

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


/*

./lshape/make_lmesh /home/kraemer/Programmieren/txts/uniform_lshape/mesh/ 2
*/









int
main(int argc, char **argv)
{
    init_h2lib(&argc, &argv);

    pstopwatch          sw;
    ptri2d              *gr_2d;
    uint                numpts, numref;
    uint                i, j;
    char                filepath[250], filename[250];
    FILE                *fptr;

    /* Parameters */
    sw = new_stopwatch();
    strcpy(filepath, argv[1]);
//    filepath = strcpy(argv[1]);
    numref = atoi(argv[2]);         /* Number of refinements*/
    assert(numref<=10);           /* 8 refs are 200000 pts*/

    (void) printf("\nCreating L-shaped mesh,\n\t%u refinement(s)\n       ", numref);

    gr_2d = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (numref + 1));
    gr_2d[0] = new_lshape_tri2d();  
//    gr_2d[0] = new_unitsquare_tri2d();  
    for (i = 0; i < numref; i++)
        gr_2d[i + 1] = refine_tri2d(gr_2d[i], NULL);
    check_tri2d(gr_2d[numref]);     
    numpts =  gr_2d[numref]->vertices;
    (void) printf("\t%u vertices\n", numpts);



 /*   (void) printf("Turning the L-shape...\n");
    for(i = 0; i < numpts; i++){
        if(gr_2d[numref]->x[i][0] > 0 && gr_2d[numref]->x[i][1] < 0){
            gr_2d[numref]->x[i][1] *= -1;
        }
    }

*/


    (void) printf("Saving mesh in file\n");
    sprintf(filename, "mesh_N%d.txt", numpts);
    strcat(filepath, filename);
    fptr = fopen(filepath, "w"); 
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







