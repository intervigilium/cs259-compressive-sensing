/******************************************************************************
 * Ray_Tracer_3D.h: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/

#ifndef __Ray_Tracer_3D__
#define __Ray_Tracer_3D__

#ifndef WIN32
#include "../config.h"
#endif

#include <math.h>
#include "define.h"
#include "Array3D.h"
#include "Vector.h"

DataType Ray_Tracer_Forward (Vector3D *source, 
                           Vector3D *detector, 
                           int N_x,
                           int N_y,
                           int N_z,
                           Array3D *image);

void Ray_Tracer_Backward (Vector3D *source,
                          Vector3D *detector,
                          int N_x,
                          int N_y,
                          int N_z,
                          DataType value,
                          Array3D *image,
                          Array3D *image_denote);

void Forward_Projection (Array3D *image,
                         DataType D,
                         int q,
                         DataType ss,
                         DataType d,
                         DataType dtheta,
                         DataType ssz,
                         int qz,
                         Array3D *sino);

void Backward_Projection (Array3D *image,
                          Array3D *image_denote,
                          DataType D,
                          int q,
                          DataType ss,
                          DataType d,
                          DataType dtheta,
                          DataType ssz,
                          int qz,
                          Array3D *sino);

void Backward_Projection2(Array3D *image,          //output 
                          Array3D *image_denote, 
                          DataType D, 
                          int q,
                          DataType ss, 
                          DataType d,
                          DataType dtheta, 
                          DataType ssz,
                          int qz, 
                          Array3D *sino);

DataType Backward_Projection_Pixel(Vector3D *source,
                                   int numberofview,
                                   DataType ss2,
                                   DataType ssz,
                                   int q,
                                   int qz,
                                   DataType ratio2,
                                   DataType D,
                                   DataType d,
                                   int N_x, 
                                   int N_y, 
                                   int N_z,
                                   int i,
                                   int j,
                                   int k,
                                   Array3D *sino);

void Backward_Projection_Column(Vector3D *source,
				   int numberofview,
				   DataType ss2,
				   DataType ssz,
				   int q,
				   int qz,
				   DataType ratio2,
				   DataType D,
				   DataType d,
                                   int N_x, 
                                   int N_y, 
                                   int N_z,
				   int i,
				   int j,
                                   Array3D *sino,
								   Array3D *image,
								   Array3D *image_denote);


void Find_Zero_Pixels (Vector3D *source,
                       Vector3D *detector,
                       Array3D  *image);

void Initialization (Array3D *image,
                     DataType D, 
                     int q, 
                     DataType ss, 
                     DataType d, 
                     DataType dtheta, 
                     DataType ssz, 
                     int qz, 
                     Array3D *sino);

#endif
