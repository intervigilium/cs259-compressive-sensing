/******************************************************************************
 * TVupdate.h: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/

#ifndef __TVupdate_H__
#define __TVupdate_H__

#include <math.h>
#include "Array3D.h"

void TVupdate(Array3D * image, Array3D * sumA, DataType alpha, int max_iter);

#endif
