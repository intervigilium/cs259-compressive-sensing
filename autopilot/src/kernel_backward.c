/*
 * Compressive Sensing kernel for FPGA implementation
 */

#include <stdlib.h>
#include <math.h>

/*
 * ************************************************
 * Function: Tracer Backward 
 * Input:  source
 *         detector
 *         image_denote 
 * Output: image  
 * Return:  
 *************************************************
 */
void Ray_Tracer_Backward(Vector3D * source,
			 Vector3D * detector,
			 int N_x,
			 int N_y,
			 int N_z,
			 DataType value,
			 Array3D * image, Array3D * image_denote)
{
	DataType L;
	Vector3D Ray_Dir;
	Vector3D Len;
	int signx, signy, signz;
	DataType absvalue_x;
	DataType absvalue_y;
	DataType absvalue_z;
	DataType lambda_min;
	DataType lambda_max;

	Vector3DInt v, u;
	int index;
	DataType x_min;
	DataType y_min;
	DataType z_min;

	DataType Tx, Ty, Tz;
	DataType length;
	int flag = 1;
	DataType lambda_0;
	Vector3D lambda;

	Vector3D_Set(&Len, 0, 0, 0);
	Vector3D_Sub(&Ray_Dir, detector, source);
	L = Vector3D_Dot(&Ray_Dir, &Ray_Dir);
	L = sqrtf(L);

	signx = (Ray_Dir.x > 0) ? 1 : -1;
	signy = (Ray_Dir.y > 0) ? 1 : -1;
	signz = (Ray_Dir.z > 0) ? 1 : -1;

	absvalue_x = ABS_VALUE(Ray_Dir.x);
	absvalue_y = ABS_VALUE(Ray_Dir.y);
	absvalue_z = ABS_VALUE(Ray_Dir.z);
	Len.x = (absvalue_x > 1.e-6) ? (L / absvalue_x) : 1.e6;
	Len.y = (absvalue_y > 1.e-6) ? (L / absvalue_y) : 1.e6;
	Len.z = (absvalue_z > 1.e-6) ? (L / absvalue_z) : 1.e6;

	lambda_min = MAX(MAX(0.0,
			     MIN(LAMBDA_X(0, source->x, detector->x, L),
				 LAMBDA_X(N_x, source->x, detector->x, L)
			     )
			 ),
			 MAX(MIN
			     (LAMBDA_Y(0, source->y, detector->y, L),
			      LAMBDA_Y(N_y, source->y, detector->y, L)),
			     MIN(LAMBDA_Z(0, source->z, detector->z, L),
				 LAMBDA_Z(N_z, source->z, detector->z, L))
			 )
	    );

	lambda_max = MIN(MIN(L,
			     MAX(LAMBDA_X(0, source->x, detector->x, L),
				 LAMBDA_X(N_x, source->x, detector->x, L)
			     )
			 ),
			 MIN(MAX
			     (LAMBDA_Y(0, source->y, detector->y, L),
			      LAMBDA_Y(N_y, source->y, detector->y, L)),
			     MAX(LAMBDA_Z(0, source->z, detector->z, L),
				 LAMBDA_Z(N_z, source->z, detector->z, L))
			 )
	    );

	if (lambda_min >= lambda_max)
		return;

	x_min = source->x + lambda_min * Ray_Dir.x / L;
	y_min = source->y + lambda_min * Ray_Dir.y / L;
	z_min = source->z + lambda_min * Ray_Dir.z / L;

	x_min = (x_min) > 0 ? x_min : 0;
	y_min = (y_min) > 0 ? y_min : 0;
	z_min = (z_min) > 0 ? z_min : 0;

	if (ABS_VALUE(x_min - (int)(x_min + 0.5)) < 1.e-6) {
		v.x = (int)(x_min + 0.5);
		u.x = v.x + signx;
		index = 1;
		if (Ray_Dir.x < 0.0)
			v.x -= 1;
	} else {
		v.x = (int)(x_min);
		if (Ray_Dir.x < 0.0)
			u.x = (int)(x_min);
		else
			u.x = (int)(x_min + 1);
	}

	if (ABS_VALUE(y_min - (int)(y_min + 0.5)) < 1.e-6) {
		v.y = (int)(y_min + 0.5);
		u.y = v.y + signy;
		index = 2;
		if (Ray_Dir.y < 0.0)
			v.y -= 1;
	} else {
		v.y = (int)(y_min);
		if (Ray_Dir.y < 0.0)
			u.y = (int)(y_min);
		else
			u.y = (int)(y_min + 1);
	}

	if (ABS_VALUE(z_min - (int)(z_min + 0.5)) < 1.e-6) {
		v.z = (int)(z_min + 0.5);
		u.z = v.z + signz;
		index = 3;
		if (Ray_Dir.z < 0.0)
			v.z -= 1;
	} else {
		v.z = (int)(z_min);
		if (Ray_Dir.z < 0.0)
			u.z = (int)(z_min);
		else
			u.z = (int)(z_min + 1);
	}

	lambda_0 = lambda_min;

	absvalue_x = ABS_VALUE(Ray_Dir.x);
	absvalue_y = ABS_VALUE(Ray_Dir.y);
	absvalue_z = ABS_VALUE(Ray_Dir.z);

	lambda.x =
	    (absvalue_x < 1e-6) ? 1e6 : LAMBDA_X(u.x, source->x, detector->x,
						 L);
	lambda.y =
	    (absvalue_y < 1e-6) ? 1e6 : LAMBDA_Y(u.y, source->y, detector->y,
						 L);
	lambda.z =
	    (absvalue_z < 1e-6) ? 1e6 : LAMBDA_Z(u.z, source->z, detector->z,
						 L);

	/// the main loop
	if (ABS_VALUE(Ray_Dir.x) > ABS_VALUE(Ray_Dir.y)) {
		if (index == 1) {
			if (MIN(lambda.z, lambda.y) > lambda_max) {
				while (1) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0)
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    Len.x * value;

					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;
				}
			}

			if (lambda.z < lambda.y) {
				Ty = lambda.y - lambda.z;
				Tz = lambda.z - lambda_0;
				Tx = lambda.x - lambda_0;
				if (Tz < Tx) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;

					}
					Tx -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}

					Tz -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;

					while (Tz >= Tx) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Tx;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Tz -= Tx;
						v.x += signx;
						if (v.x >= N_x || v.x < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}

					Tx -= Tz;
					Tz = Len.z;

					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				}
			} else {
				Tz = lambda.z - lambda.y;
				Ty = lambda.y - lambda_0;
				Tx = lambda.x - lambda_0;
				if (Ty < Tx) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;

					while (Ty >= Tx) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Tx;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Ty -= Tx;
						v.x += signx;
						if (v.x >= N_x || v.x < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;
				}
			}
		} else if (index == 2) {
			Ty = Len.y;
			Tz = lambda.z - lambda_0;
			Tx = lambda.x - lambda_0;
		} else {
			Tz = Len.z;
			Ty = lambda.y - lambda_0;
			Tx = lambda.x - lambda_0;
		}

		while (flag) {
			if (Tz < Ty) {
				Ty -= Tz;
				if (Tz < Tx) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tz -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;

					while (Tz >= Tx) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Tx;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}

						Tz -= Tx;
						v.x += signx;
						if (v.x >= N_x || v.x < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				}
			} else {
				Tz -= Ty;
				if (Ty < Tx) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;

					while (Ty >= Tx) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Tx;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Ty -= Tx;
						v.x += signx;
						if (v.x >= N_x || v.x < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;
				}
			}
		}
	} else {
		if (index == 2) {
			if (MIN(lambda.z, lambda.x) > lambda_max) {
				while (1) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0)
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    Len.y * value;

					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;
				}
			}
			if (lambda.z < lambda.x) {
				Tx = lambda.x - lambda.z;
				Tz = lambda.z - lambda_0;
				Ty = lambda.y - lambda_0;
				if (Tz < Ty) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tz -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;

					while (Tz >= Ty) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Ty;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Tz -= Ty;
						v.y += signy;
						if (v.y >= N_y || v.y < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				}
			} else {
				Tz = lambda.z - lambda.x;
				Tx = lambda.x - lambda_0;
				Ty = lambda.y - lambda_0;
				if (Tx < Ty) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;

					while (Tx >= Ty) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Ty;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Tx -= Ty;
						v.y += signy;
						if (v.y >= N_y || v.y < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;
				}
			}
		} else if (index == 1) {
			Tx = Len.x;
			Tz = lambda.z - lambda_0;
			Ty = lambda.y - lambda_0;
		} else {
			Tz = Len.z;
			Tx = lambda.x - lambda_0;
			Ty = lambda.y - lambda_0;
		}

		while (1) {
			if (Tz < Tx) {
				Tx -= Tz;
				if (Tz < Ty) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tz -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;

					while (Tz >= Ty) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Ty;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Tz -= Ty;
						v.y += signy;
						if (v.y >= N_y || v.y < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tz;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tz;
					Tz = Len.z;
					v.z += signz;
					if (v.z >= N_z || v.z < 0)
						return;
				}
			} else {
				Tz -= Tx;
				if (Tx < Ty) {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty = Ty - Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;
				} else {
					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Ty;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Tx -= Ty;
					Ty = Len.y;
					v.y += signy;
					if (v.y >= N_y || v.y < 0)
						return;

					while (Tx >= Ty) {
						if (ARRAY3DData
						    (image_denote, v.x, v.y,
						     v.z) > 0) {
							length = Ty;
							ARRAY3DData(image, v.x,
								    v.y, v.z) +=
							    length * value;
						}
						Tx -= Ty;
						v.y += signy;
						if (v.y >= N_y || v.y < 0)
							return;
					}

					if (ARRAY3DData
					    (image_denote, v.x, v.y, v.z) > 0) {
						length = Tx;
						ARRAY3DData(image, v.x, v.y,
							    v.z) +=
						    length * value;
					}
					Ty -= Tx;
					Tx = Len.x;
					v.x += signx;
					if (v.x >= N_x || v.x < 0)
						return;
				}
			}
		}
	}
}

void Backward_Projection(Array3D * image,	//output 
			 Array3D * image_denote,
			 DataType D,
			 int q,
			 DataType ss,
			 DataType d,
			 DataType dtheta, DataType ssz, int qz, Array3D * sino)
{
	int N_x = image->index1Size;
	int N_y = image->index2Size;
	int N_z = image->index3Size;
	int NS = (int)(359.0 / dtheta);
	DataType thetaS = 0.;
	DataType dtheta2 = (DataType) (dtheta * PI / 180.);
	DataType ss2 = (DataType) (ss * PI / 180.);
	int i, j, k;
	Vector3D source, detector;

	Array3D_initialize_value(image->dataPtr, 0.0, N_x, N_y, N_z);

	for (i = 0; i <= NS; i++) {
		source.x = D * cos(thetaS) + d;
		source.y = D * sin(thetaS) + d;
		source.z = d;
		for (j = -q; j < q + 1; j++) {
			detector.x = D * cos(j * ss2 + PI + thetaS) + d;
			detector.y = D * sin(j * ss2 + PI + thetaS) + d;
			for (k = -qz; k < qz + 1; k++) {
				DataType value =
				    ARRAY3DData(sino, (j + q), (k + qz), i);
				detector.z = k * ssz + d;

				if (value > 1.e-6)
					Ray_Tracer_Backward(&source, &detector,
							    N_x, N_y, N_z,
							    value, image,
							    image_denote);
			}
		}
		thetaS += dtheta2;
	}
}
