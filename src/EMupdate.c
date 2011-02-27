/******************************************************************************
 * EMupdate.cpp: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/

#include "EMupdate.h"

void EMupdate(Array3D * sino,
	      Array3D * image,
	      Array3D * image_denote,
	      Array3D * sumA,
	      int max_iter,
	      DataType D,
	      int q,
	      DataType ss,
	      DataType d, DataType dtheta, DataType ssz, int qz)
{
    int M_x = sino->index1Size;
    int M_y = sino->index2Size;
    int M_z = sino->index3Size;
    int N_x = image->index1Size;
    int N_y = image->index2Size;
    int N_z = image->index3Size;
    int i, j, k, l;

    Array3D Af;
    Array3D gAf;
    Array3D AtAf;

    Array3D_malloc(&Af, M_x, M_y, M_z);
    Array3D_malloc(&gAf, M_x, M_y, M_z);
    Array3D_malloc(&AtAf, N_x, N_y, N_z);

    for (i = 0; i < max_iter; i++) {
	Forward_Projection(image, D, q, ss, d, dtheta, ssz, qz, &Af);

	for (j = 0; j < M_x; j++)
	    for (k = 0; k < M_y; k++)
		for (l = 0; l < M_z; l++)
		    if (ARRAY3DData((&Af), j, k, l) > 1.e-6) {
			ARRAY3DData((&gAf), j, k, l) =
			    ARRAY3DData(sino, j, k, l) / ARRAY3DData((&Af),
								     j, k,
								     l);
		    }

	DBG1(fprintf(stderr, "\n");
	    )

	    Backward_Projection(&AtAf, image_denote, D, q, ss, d, dtheta,
				ssz, qz, &gAf);

	DBG1(fprintf(stderr, "\n");
	    )

	    for (j = 0; j < N_x; j++)
	    for (k = 0; k < N_y; k++)
		for (l = 0; l < N_z; l++) {
		    ARRAY3DData(image, j, k, l) =
			ARRAY3DData(image, j, k, l) * ARRAY3DData((&AtAf),
								  j, k,
								  l) /
			ARRAY3DData(sumA, j, k, l);
		}
    }

    Array3D_free(&Af);
    Array3D_free(&gAf);
    Array3D_free(&AtAf);
}
