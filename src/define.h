#ifndef __DEFINE_H__
#define __DEFINE_H__

#define DEBUG_CHECK
#define DEBUG_FILEOUTPUT

#define DEBUG_LEVEL1
#ifdef  DEBUG_LEVEL1
#define DBG1(expr)  expr
#else
#define DBG1(expr)
#endif

//#define DEBUG_LEVEL2
#ifdef  DEBUG_LEVEL2
#define DBG2(expr)  expr
#else
#define DBG2(expr)
#endif

#define SINGLEPOINT
#ifdef SINGLEPOINT
#define DataType float
#define SQRT  sqrtf
#define FLOOR floorf
#else
#define DataType double
#define SQRT  sqrt
#define FLOOR floor
#endif

#define PI (4.*atan2(1.0,1.0))

#ifndef MAX
#define MAX( x, y ) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef MIN
#define MIN( x, y ) ( ((x) < (y)) ? (x) : (y) )
#endif

#ifndef MAX3
#define MAX3(x,y,z) MAX(MAX(x,y),z)
#endif

#ifndef MIN3
#define MIN3(x,y,z) MIN(MIN(x,y),z)
#endif

#define ABS_VALUE(x) ( (x < 0) ? -(x) : (x) )

#define LAMBDA_X(i, x_s, x_d, L) (L*((DataType)i-x_s)/(x_d-x_s))
#define LAMBDA_Y(j, y_s, y_d, L) (L*((DataType)j-y_s)/(y_d-y_s))
#define LAMBDA_Z(k, z_s, z_d, L) (L*((DataType)k-z_s)/(z_d-z_s))


#endif
