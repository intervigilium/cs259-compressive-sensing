/******************************************************************************
 * TVupdate.cpp: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/

#include <math.h>

#include "define.h"
#include "TVupdate.h"

void TVupdate(Array3D * image,
	      Array3D * sumA, DataType alpha, int max_iter)
{
    int N_x = image->index1Size - 1;
    int N_y = image->index2Size - 1;
    int N_z = image->index3Size - 1;

    DataType c1, c2, c3, c4;
    DataType eps = (DataType) (1e-6);
    int i, j, k, l;

    for (i = 0; i < max_iter; i++) {
	for (j = 1; j < N_x; j++)
	    for (k = 1; k < N_y; k++)
		for (l = 1; l < N_z; l++) {
		    DataType data000;
		    DataType data100, data010, data001;
		    DataType datan100, datan110, datan101;
		    DataType data1n10, data0n10, data0n11;
		    DataType data10n1, data00n1, data01n1;
		    DataType sum;
		    DataType sub100, sub010, sub001;
		    DataType sub000, subn110, subn101;
		    DataType sub1n10, sub0002, sub0n11;
		    DataType sub10n1, sub01n1, sub0003;
		    DataType temp;

		    data000 = ARRAY3DData(image, j, k, l);
		    if (data000 > 1e-6) {
			data100 = ARRAY3DData(image, (j + 1), k, l);
			data010 = ARRAY3DData(image, j, (k + 1), l);
			data001 = ARRAY3DData(image, j, k, (l + 1));
			datan100 = ARRAY3DData(image, (j - 1), k, l);
			datan110 = ARRAY3DData(image, (j - 1), (k + 1), l);
			datan101 = ARRAY3DData(image, (j - 1), k, (l + 1));
			data1n10 = ARRAY3DData(image, (j + 1), (k - 1), l);
			data0n10 = ARRAY3DData(image, j, (k - 1), l);
			data0n11 = ARRAY3DData(image, j, (k - 1), (l + 1));
			data10n1 = ARRAY3DData(image, (j + 1), k, (l - 1));
			data00n1 = ARRAY3DData(image, j, k, (l - 1));
			data01n1 = ARRAY3DData(image, j, (k + 1), (l - 1));
			sum = ARRAY3DData(sumA, j, k, l);


			sub100 = data100 - data000;
			sub010 = data010 - data000;
			sub001 = data001 - data000;


			sub000 = data000 - datan100;
			subn110 = datan110 - datan100;
			subn101 = datan101 - datan100;


			sub1n10 = data1n10 - data0n10;
			sub0002 = data000 - data0n10;
			sub0n11 = data0n11 - data0n10;

			sub10n1 = data10n1 - data00n1;
			sub01n1 = data01n1 - data00n1;
			sub0003 = data000 - data00n1;

			c1 = data000 / SQRT(eps + sub100 * sub100 +
					    sub010 * sub010 +
					    sub001 * sub001) / sum;
			c2 = data000 / SQRT(eps + sub000 * sub000 +
					    subn110 * subn110 +
					    subn101 * subn101) / sum;
			c3 = data000 / SQRT(eps + sub1n10 * sub1n10 +
					    sub0002 * sub0002 +
					    sub0n11 * sub0n11) / sum;
			c4 = data000 / SQRT(eps + sub10n1 * sub10n1 +
					    sub01n1 * sub01n1 +
					    sub0003 * sub0003) / sum;


			temp = alpha * ARRAY3DData(image, j, k, l)
			    + c1 * (ARRAY3DData(image, (j + 1), k, l) +
				    ARRAY3DData(image, j, (k + 1),
						l) + ARRAY3DData(image, j,
								 k,
								 (l + 1)))
			    + c2 * ARRAY3DData(image, (j - 1), k, l)
			    + c3 * ARRAY3DData(image, j, (k - 1), l)
			    + c4 * ARRAY3DData(image, j, k, (l - 1));

			ARRAY3DData(image, j, k, l) =
			    (DataType) (temp /
					(alpha + 3.0 * c1 + c2 + c3 + c4));
		    }
/*
if (image(j,k,l) > 1e-8)
{
c1 =   image(j,k,l)/sqrt(eps + pow(image(j + 1,k,l) - image(j,k,l),2) + pow(image(j,k + 1,l) - image(j,k,l),2) + pow(image(j,k,l + 1) - image(j,k,l),2))/sumA(j,k,l);
c2 = image(j,k,l)/sqrt(eps + pow(image(j,k,l) - image(j - 1,k,l),2) + pow(image(j - 1,k + 1,l) - image(j - 1,k,l),2) + pow(image(j - 1,k,l + 1) - image(j - 1,k,l),2))/sumA(j,k,l);
c3 = image(j,k,l)/sqrt(eps + pow(image(j + 1,k - 1,l) - image(j,k - 1,l),2) + pow(image(j,k,l) - image(j,k - 1,l),2) + pow(image(j,k - 1,l + 1) - image(j,k - 1,l),2))/sumA(j,k,l);
c4 = image(j,k,l)/sqrt(eps + pow(image(j + 1,k,l - 1) - image(j,k,l - 1),2) + pow(image(j,k + 1,l - 1) - image(j,k,l - 1),2) + pow(image(j,k,l) - image(j,k,l - 1),2))/sumA(j,k,l);
image(j,k,l) = (alpha*image(j,k,l) + c1*(image(j + 1,k,l) + image(j,k + 1,l) + image(j,k,l + 1)) + c2*image(j - 1,k,l) + c3*image(j,k - 1,l) + c4*image(j,k,l - 1))/(alpha + 3.0*c1 + c2 + c3 + c4);				
}				
*/
		};

	for (j = 1; j < N_x; j++)
	    for (k = 1; k < N_y; k++) {
		ARRAY3DData(image, j, k, 0) = ARRAY3DData(image, j, k, 1);
		ARRAY3DData(image, j, k, N_z) =
		    ARRAY3DData(image, j, k, (N_z - 1));
	    }
    }
}
