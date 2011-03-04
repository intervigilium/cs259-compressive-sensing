/******************************************************************************
 * EMupdate.h: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/

#ifndef __EMupdate_H__
#define __EMupdate_H__

#include "Array3D.h"
#include "Ray_Tracer_3D.h"

void EMupdate(Array3D * sino,
	      Array3D * image,
	      Array3D * image_denote,
	      Array3D * sumA,
	      int max_iter,
	      DataType D,
	      int q,
	      DataType ss, DataType d, DataType dtheta, DataType ssz, int qz);

#endif
