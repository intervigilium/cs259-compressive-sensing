#include "define.h"
#include "Array3D.h"
#include "Vector.h"
#include "Ray_Tracer_3D.h"
#include "EMupdate.h"
#include "TVupdate.h"
#include "EMTV3D.h"


int CPU_Routine(Array3D * image, Array3D * sumA, Array3D * image_denote,
		Array3D * sino)
{
    int N_x = 128, N_y = 128, N_z = 128;	// the size of the image, N_x * N_y * N_z
    DataType D = 125, d = 64;	// the parameters used to define the geometry of the machine
    DataType theta = 10.0;
    int q = 150;
    int qz = 128;
    DataType ss = 0.5;
    DataType ssz = 1.0;
    int max_iter_EM = 3;
    int max_iter_TV = 10;
    DataType alpha = 5;
    int iter;
    int64_t i_start_whole, i_end_whole;
    int64_t i_start, i_end;

#ifdef DEBUG_FILEOUTPUT
    char filename[100];
#endif

#ifdef DEBUG_FILEOUTPUT
    writeImgArray3D("../data/cpu_readinput.dat", image);
#endif

    i_start_whole = emtv3d_timer_s();

    fprintf(stderr, "CPU Backward_Projection begins \n");
    i_start = emtv3d_timer_us();
    Backward_Projection(sumA, image_denote, D, q, ss, d, theta, ssz, qz, sino);	// find summation of columns
    i_end = emtv3d_timer_us();
    fprintf(stderr, "\nCPU Backward_Projection uses %lld us\n",
	    (i_end - i_start));

#ifdef DEBUG_FILEOUTPUT
    writeImgArray3D("../data/cpu_sumA.dat", sumA);
#endif

#if 0
#endif

    fprintf(stderr, "CPU Forward_Projection begins \n");
    i_start = emtv3d_timer_us();
    Forward_Projection(image, D, q, ss, d, theta, ssz, qz, sino);	// creat the sinogram data, for final use, load the sinogram data
    i_end = emtv3d_timer_us();
    fprintf(stderr, "\n CPU Forward_Projection uses %lld us\n",
	    (i_end - i_start));

#ifdef DEBUG_FILEOUTPUT
    writeImgArray3D("../data/cpu_sino.dat", sino);
    //writeImgArray3DWithIndex("../data/sinowithindex.dat",sino);
#endif

    fprintf(stderr, "CPU Initialization\n");
    Initialization(image_denote, D, q, ss, d, theta, ssz, qz, sino);	// initialize the image, make sure some pixels are zeros;
    Array3D_copy(image->dataPtr, image_denote->dataPtr, N_x, N_y, N_z);

#ifdef DEBUG_FILEOUTPUT
    writeImgArray3D("../data/cpu_initialization.dat", image);
#endif

    fprintf(stderr, "\n");
    i_start = emtv3d_timer_s();
    for (iter = 0; iter < 100; iter++) {
	fprintf(stderr, "CPU EMupdate %d \n", iter);
	EMupdate(sino, image, image_denote, sumA, max_iter_EM, D, q, ss, d, theta, ssz, qz);	// EMupdating

#ifdef DEBUG_FILEOUTPUT
	sprintf(filename, "../data/cpu_EMupdate%d.dat", iter);
	writeImgArray3D(filename, image);
#endif
	fprintf(stderr, "CPU TVupdate %d \n", iter);
	TVupdate(image, sumA, alpha, max_iter_TV);

#ifdef DEBUG_FILEOUTPUT
	sprintf(filename, "../data/cpu_TVupdate%d.dat", iter);
	writeImgArray3D(filename, image);
#endif
    }
    i_end = emtv3d_timer_s();
    i_end_whole = emtv3d_timer_s();

    fprintf(stderr, "CPU \nEMTV update uses %lld seconds\n",
	    (i_end - i_start));
    fprintf(stderr, "CPU total using %lld seconds\n",
	    (i_end_whole - i_start_whole));

    return 0;
}
