/******************************************************************************
 * EMTV3D.cpp: EMTV for 3D image reconstruction
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
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>


#ifdef WIN32
#include "extras/getopt.h"
#else
#include <sys/time.h>
#include <getopt.h>
#endif

#ifndef WIN32
#include "../config.h"
#endif


#include "define.h"
#include "Array3D.h"
#include "Vector.h"
#include "Ray_Tracer_3D.h"
#include "EMupdate.h"
#include "TVupdate.h"
#include "CPU_Routine.h"
#include "EMTV3D.h"
#include "papi.h"


char inputfilename[30] = "../data/phantom_data128";	//default input  file name
char outputfilename[30] = "../data/rec_image_data";	//default output file name
char outputfilename_gpu[30] = "../data/rec_image_data_gpu";	//default output file name
unsigned batch_id = 0;

int64_t emtv3d_timer_us(void)
{
#ifdef WIN32
    struct _timeb tb;
    _ftime(&tb);
    return ((int64_t) tb.time * (1000) + (int64_t) tb.millitm) * (1000);
#else
    struct timeval tv_date;
    gettimeofday(&tv_date, NULL);
    return ((int64_t) tv_date.tv_sec * 1000000 +
	    (int64_t) tv_date.tv_usec);
#endif
}

int64_t emtv3d_timer_s(void)
{
#ifdef WIN32
    struct _timeb tb;
    _ftime(&tb);
    return ((int64_t) tb.time);
#else
    struct timeval tv_date;
    gettimeofday(&tv_date, NULL);
    return ((int64_t) tv_date.tv_sec);
#endif
}


void readImgArray3D(char *in_file, Array3D * image)
{
    int N_x = image->index1Size;
    int N_y = image->index2Size;
    int N_z = image->index3Size;
    FILE *inputfile;
    int i, j, k;

    inputfile = fopen(in_file, "r");
    if (inputfile == NULL) {
	printf("Cannot open the input file!\n");
	exit(1);
    }

    for (i = 0; i < N_x; i++)
	for (j = 0; j < N_y; j++)
	    for (k = 0; k < N_z; k++) {
#ifdef SINGLEPOINT
		fscanf(inputfile, "%f", ARRAY3DAddr(image, i, j, k));
#else
		fscanf(inputfile, "%lf", ARRAY3DAddr(image, i, j, k));
#endif
	    }
    fclose(inputfile);
}

void writeImgArray3D(char *out_file, Array3D * image)
{
    int N_x = image->index1Size;
    int N_y = image->index2Size;
    int N_z = image->index3Size;
    FILE *outputfile;
    int i, j, k;

    outputfile = fopen(out_file, "w");
    if (outputfile == NULL) {
	printf("Cannot open the output file!\n");
	exit(1);
    }
    for (i = 0; i < N_x; i++)
	for (j = 0; j < N_y; j++)
	    for (k = 0; k < N_z; k++) {
#ifdef SINGLEPOINT
		fprintf(outputfile, "%f\n", ARRAY3DData(image, i, j, k));
#else
		fprintf(outputfile, "%lf\n", ARRAY3DData(image, i, j, k));
#endif
	    };
    fclose(outputfile);
}

void writeImgArray3DWithIndex(char *out_file, Array3D * image)
{
    int N_x = image->index1Size;
    int N_y = image->index2Size;
    int N_z = image->index3Size;
    FILE *outputfile;
    int i, j, k;

    outputfile = fopen(out_file, "w");
    if (outputfile == NULL) {
	printf("Cannot open the output file!\n");
	exit(1);
    }
    for (i = 0; i < N_x; i++)
	for (j = 0; j < N_y; j++)
	    for (k = 0; k < N_z; k++) {
#ifdef SINGLEPOINT
		fprintf(outputfile, "%3d, %3d, %3d, %f\n", i, j, k,
			ARRAY3DData(image, i, j, k));
#else
		fprintf(outputfile, "%3d, %3d, %3d, %lf\n", i, j, k,
			ARRAY3DData(image, i, j, k));
#endif
	    };
    fclose(outputfile);
}

void CompareVector(DataType * data0, DataType * data1, int size)
{
    int i;
    DataType diff = 0.0;
    DataType diffall = 0.0;

    for (i = 0; i < size; i++) {
	diff = data0[i] - data1[i];

	if (diff > 0.1)
#ifdef SINGLEPOINT
	    fprintf(stderr, "%d, data0[%d] = %f data1[%d]= %f \n", i, i,
		    data0[i], i, data1[i]);
#else
	    fprintf(stderr, "%d, data0[%d] = %lf data1[%d]= %lf \n", i, i,
		    data0[i], i, data1[i]);
#endif

	diffall += diff;
    }

#ifdef SINGLEPOINT
    fprintf(stderr, "diffall = %f\n", diffall);
#else
    fprintf(stderr, "diffall = %lf\n", diffall);
#endif
}

int help()
{
    printf("./emtv3d -i inputfilename -o outputfilename\n");
    return 0;
}


int configure(int argc, char *argv[])
{
    char option_def[] = "?vh:i:o:b:";
    int ch;

    for (;;) {
	ch = getopt(argc, argv, option_def);
	if (ch == -1)
	    break;

	switch (ch) {
	case '?':
	case 'v':
	case 'h':
	    help();
	    return -1;
	case 'i':
	    strcpy(inputfilename, optarg);
	    printf("input filename : %s \n", inputfilename);
	    break;
	case 'o':
	    strcpy(outputfilename, optarg);
	    printf("output filename : %s \n", outputfilename);
	    break;
	case 'b':
	    batch_id = atoi(optarg);
	    printf("papi batch id is : %d \n", batch_id);
	    break;
	default:
	    printf("Unknow option detection\n");
	    return 0;
	}
    }

    return 0;
}


int main(int argc, char **argv)
{
    int N_x = 128, N_y = 128, N_z = 128;	// the size of the image, N_x * N_y * N_z
    DataType theta = 10.0;
    int q = 150;
    int qz = 128;

    Array3D cpu_image;		// used to store the reconstructed image
    Array3D cpu_sumA;		// used to store the summation of columns
    Array3D cpu_image_denote;
    Array3D cpu_sino;		// used to store the sinogram data

    if (configure(argc, argv) == -1)
	return -1;

#ifdef SINGLEPOINT
    fprintf(stdout, "Using single point data-type\n");
#else
    fprintf(stdout, "Using double point data-type\n");
#endif
    fprintf(stdout, "input image: %s\n", inputfilename);
    fprintf(stdout, "output image: %s\n", outputfilename);

    //initial global buffers
    Array3D_malloc(&cpu_image, N_x, N_y, N_z);
    Array3D_malloc(&cpu_sumA, N_x, N_y, N_z);
    Array3D_malloc(&cpu_image_denote, N_x, N_y, N_z);
    Array3D_initialize_value((cpu_image_denote.dataPtr), 1.0, N_x, N_y,
			     N_z);
    Array3D_malloc(&cpu_sino, 2 * q + 1, 2 * qz + 1,
		   (int) floorf(359. / theta) + 1);
    Array3D_initialize_value((cpu_sino.dataPtr), 1.0, 2 * q + 1,
			     2 * qz + 1, (int) floorf(359. / theta) + 1);

    //load input data 
    readImgArray3D(inputfilename, &(cpu_image));
    //CPU routine 

    int Events[5];
    u_long_long papi_values[5];
    util_start_papi(batch_id, Events);
    CPU_Routine(&cpu_image, &cpu_sumA, &cpu_image_denote, &cpu_sino);
    util_stop_papi(batch_id, papi_values);
    util_print_papi(batch_id, papi_values, (batch_id == 0));
    //output CPU results
    writeImgArray3D(outputfilename, &cpu_image);

    //free the global buffers
    Array3D_free(&cpu_image);
    Array3D_free(&cpu_sumA);
    Array3D_free(&cpu_image_denote);
    Array3D_free(&cpu_sino);

    return 0;
}
