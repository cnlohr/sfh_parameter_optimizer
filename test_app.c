#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sfh_parameter_optimizer.h"

typedef struct
{
	double uvTarget[2];
	double uvImage[2];
} cal_point;

typedef struct
{
	int nrcalpts;
	cal_point * calpts;

	int image_width;
	int image_height;
	int image_major_axis;

	// Camera Extrinsics
	// Try euler angles
	double rotation[3];
	double position[3];

	double m34extrinsics[12]; // Transform from view to world.

	double thisImageError;
} cal_image;

int xcodes, ycodes, axismaxcodes;

#include "debug_lenscal_opt.h"

#define CAM_CAL_ORDER_K 6
#define CAM_CAL_ORDER_P 2
#define CAM_CAL_ORDER_S 4

// https://en.wikipedia.org/wiki/Distortion_(optics)
// 
typedef struct
{
	double Kn[CAM_CAL_ORDER_K];
	double Pn[CAM_CAL_ORDER_P];
	double Sn[CAM_CAL_ORDER_S];

	double cx, cy; // Distortion center.
	double fx, fy;

	double m34intrinsics[12];  // Transform from image space to view space.

} camcal_t;

#define OPT_RMS

#ifdef OPT_RMS
#define SERR( x ) (sqrt(x))
#else
#define SERR( x ) ( x )
#endif

camcal_t camcal;

#define m00 0
#define m01 1
#define m02 2
#define m03 3
#define m10 4
#define m11 5
#define m12 6
#define m13 7
#define m20 8
#define m21 9
#define m22 10
#define m23 11
#define m30 12
#define m31 13
#define m32 14
#define m33 15

void GenerateRotationMatrix( double * outmat, double * hpr )
{
#if 0
	// https://mathworld.wolfram.com/EulerAngles.html
	// Doesn't work well (first option from above site)
	float sA = sin( hpr[0] );
	float cA = cos( hpr[0] );
	float sB = sin( hpr[1] );
	float cB = cos( hpr[1] );
	float sC = sin( hpr[2] );
	float cC = cos( hpr[2] );
	outmat[m00] =   cC*cA - cB*sA*sC;
	outmat[m01] =   cC*sA + cB*cA*sC;
	outmat[m02] =   sC*sB;
	outmat[m03] = 0;
	outmat[m10] = - sC*cA - cB*sA*cC;
	outmat[m11] = - sC*sA + cB*cA*cC;
	outmat[m12] =   cC*sB;
	outmat[m13] = 0;
	outmat[m20] =   sB*sA;
	outmat[m21] = - sB*cA;
	outmat[m22] =   cB;
	outmat[m23] = 0;
#else
	// https://mathworld.wolfram.com/EulerAngles.html
	float sA = sin( hpr[0] );
	float cA = cos( hpr[0] );
	float sB = sin( hpr[1] );
	float cB = cos( hpr[1] );
	float sC = sin( hpr[2] );
	float cC = cos( hpr[2] );
	outmat[m00] = cB*cA;
	outmat[m01] = cB*sA;
	outmat[m02] = -sB;
	outmat[m03] = 0;
	outmat[m10] = sC*sB*cA-cC*sA;
	outmat[m11] = sC*sB*sA+cC*cA;
	outmat[m12] = cB*sC;
	outmat[m13] = 0;
	outmat[m20] = cC*sB*cA+sC*sA;
	outmat[m21] = cC*sB*sA-sC*cA;
	outmat[m22] = cB*cC;
	outmat[m23] = 0;
#endif
}


// Mostly from cnovrmath
void matrix34ptransform( double * pout, const double * pin, const double * f )
{
	double ptmp[2];
	ptmp[0] = pin[0] * f[m00] + pin[1] * f[m01] + pin[2] * f[m02] + f[m03];
	ptmp[1] = pin[0] * f[m10] + pin[1] * f[m11] + pin[2] * f[m12] + f[m13];
	pout[2] = pin[0] * f[m20] + pin[1] * f[m21] + pin[2] * f[m22] + f[m23];
	pout[0] = ptmp[0];
	pout[1] = ptmp[1];
}

double square( double s )
{
	return s*s;
}

double terror( double d )
{
	if( d < 0 ) return -d;
	else return d;
}

void ComputeErrorLocal( camcal_t * cc, cal_image * image )
{
	int i;
	double totalError = 0.0;

	GenerateRotationMatrix( image->m34extrinsics, image->rotation );
	double * ip = image->position;
	image->m34extrinsics[m03] = ip[0];
	image->m34extrinsics[m13] = ip[1];
	image->m34extrinsics[m23] = ip[2];

	cc->m34intrinsics[m00] = cc->fx;
	cc->m34intrinsics[m01] = 0;
	cc->m34intrinsics[m02] = cc->cx;
	cc->m34intrinsics[m03] = 0;
	cc->m34intrinsics[m10] = 0;
	cc->m34intrinsics[m11] = cc->fy;
	cc->m34intrinsics[m12] = cc->cy;
	cc->m34intrinsics[m13] = 0;
	cc->m34intrinsics[m20] = 0;
	cc->m34intrinsics[m21] = 0;
	cc->m34intrinsics[m22] = 1;
	cc->m34intrinsics[m23] = 0;

	for( i = 0; i < image->nrcalpts; i++ )
	{
		cal_point * pt = image->calpts + i;

		// Next, apply rotation matrix
		double pOut[3];

		// XXX TODO: this should be , 1. no? Otherwise center is not gathered.
		double p[3] = { -pt->uvTarget[0], -pt->uvTarget[1], -1 };

		matrix34ptransform( pOut, p, image->m34extrinsics );
		matrix34ptransform( p, pOut, cc->m34intrinsics );

		// TODO: Is this divide-by-z ok?
		double pCalc[2] = { p[0]/p[2], p[1]/p[2] };

		float rsq = pCalc[0]*pCalc[0] + pCalc[1]*pCalc[1];

		double pKterm = 
			(1+cc->Kn[0]*rsq+cc->Kn[1]*rsq*rsq+cc->Kn[2]*rsq*rsq*rsq)/
			(1+cc->Kn[3]*rsq+cc->Kn[4]*rsq*rsq+cc->Kn[5]*rsq*rsq*rsq);

		double pOx = pCalc[0]*pKterm + 2*cc->Pn[0]*pCalc[0]*pCalc[1] + cc->Pn[1]*(rsq+2*pCalc[0]*pCalc[0]) + cc->Sn[0]*rsq+cc->Sn[1]*rsq*rsq;
		double pOy = pCalc[1]*pKterm + 2*cc->Pn[1]*pCalc[0]*pCalc[1] + cc->Pn[0]*(rsq+2*pCalc[1]*pCalc[1]) + cc->Sn[2]*rsq+cc->Sn[3]*rsq*rsq;

		pCalc[0] = pOx;
		pCalc[1] = pOy;
		double targetX = pt->uvImage[0];
		double targetY = pt->uvImage[1];

#ifdef OPT_RMS
		double err = square(pCalc[0]-targetX) + square(pCalc[1]-targetY);
#else
		double err = terror(pCalc[0]-targetX) + terror(pCalc[1]-targetY);
#endif
//exit( 0 );
		totalError += err;
	}

	image->thisImageError = totalError;
}

void SetCalibrationDefaults( camcal_t * cc, cal_image * calimages, int nrcalimages )
{
	int i;
	for( i = 0; i < nrcalimages; i++ )
	{
		cal_image * ci = calimages + i;
		ci->rotation[0] = 0;
		ci->rotation[1] = 0;
		ci->rotation[2] = 0;
		ci->position[0] = 0;
		ci->position[1] = 0;
		ci->position[2] =-1;
	}

	// HMM Is this the inverse?!
	cc->cx = 0.0;
	cc->cy = 0.0;
	cc->fx = 1;
	cc->fy = 1;
}

void AdjustLocal( void * vOpaque, int iOpaque, double v )
{
	*((double*)vOpaque) = v;
	ComputeErrorLocal( &camcal, calimages + iOpaque );
}

void AdjustGlobal( void * vOpaque, int iOpaque, double v )
{
	*((double*)vOpaque) = v;
	int i;
	for( i = 0; i < nrcalimages; i++ )
	{
		ComputeErrorLocal( &camcal, calimages + i );
	}
}

double GetGlobalError( void * vOpaque, int iOpaque )
{
	double totalError = 0;
	double totalPoints = 0;
	int i;
	for( i = 0; i < nrcalimages; i++ )
	{
		cal_image * c = calimages + i;
		totalError += c->thisImageError;
		totalPoints += c->nrcalpts;
	}

	return SERR( totalError / totalPoints );
}

int ContinueGoing( struct optimum_ctx_t * ctx )
{
	if( ( ctx->iteration % 100 ) == 0 )
		printf( "%d %f\n", ctx->iteration, ctx->currentError );
	return ctx->iteration < 20000;
}

double  Calibrate( camcal_t * cc, cal_image * calimages, int nrcalimages )
{
	// Also read https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/tr98-71.pdf
	// Notes from https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
	//
	// Step 1: Solve for K and p.
	//
	// (u, v) = K * p
	// Which means:
	//
	//   [ u ]   [ fx  0 cx ] [ r11 r12 r13 tx ] [ Xw ]
	// s [ v ] = [  0 fy cy ] [ r21 r22 r23 ty ] [ Yw ]
	//   [ 1 ]   [  0  0  1 ] [ r31 r32 r33 tz ] [ Zw ]
	//                                           [  1 ]
	// Step 2: Ease into solving for Kn, Kp.
	//
	// Step 3: Eventually also solve for the very fine lens correction parameters.

	SetCalibrationDefaults( cc, calimages, nrcalimages );

	int tuningVars = 16 + 6 * nrcalimages;

	optimum_variable variables[tuningVars];

	int i;
	for( i = 0; i < nrcalimages; i++ )
	{
		cal_image * ci = calimages + i;
		variables[16+i*6+0] = (optimum_variable){ AdjustLocal, &ci->rotation[0], i, 0.5, 0.0, 0.2 };
		variables[16+i*6+1] = (optimum_variable){ AdjustLocal, &ci->rotation[1], i, 0.5, 0.0, 0.2 };
		variables[16+i*6+2] = (optimum_variable){ AdjustLocal, &ci->rotation[2], i, 0.5, 0.0, 0.2 };
		variables[16+i*6+3] = (optimum_variable){ AdjustLocal, &ci->position[0], i, 0.5, 0.0, 0.2 };
		variables[16+i*6+4] = (optimum_variable){ AdjustLocal, &ci->position[1], i, 0.5, 0.0, 0.2 };
		variables[16+i*6+5] = (optimum_variable){ AdjustLocal, &ci->position[2], i, 0.5, 0.0, 0.2 };
	}

	variables[ 0] = (optimum_variable){ AdjustGlobal, &camcal.fx, 0, 0.04, 1.0, 0.001 };
	variables[ 1] = (optimum_variable){ AdjustGlobal, &camcal.fy, 0, 0.04, 1.0, 0.001 };
	variables[ 2] = (optimum_variable){ AdjustGlobal, &camcal.cx, 0, 0.01, 0.0, 0.001 };
	variables[ 3] = (optimum_variable){ AdjustGlobal, &camcal.cy, 0, 0.01, 0.0, 0.001 };

	variables[ 4] = (optimum_variable){ AdjustGlobal, &camcal.Kn[0], 0, 0.01, 0.0, 0.001 };
	variables[ 5] = (optimum_variable){ AdjustGlobal, &camcal.Kn[1], 0, 0.01, 0.0, 0.001 };
	variables[ 6] = (optimum_variable){ AdjustGlobal, &camcal.Kn[2], 0, 0.01, 0.0, 0.001 };
	variables[ 7] = (optimum_variable){ AdjustGlobal, &camcal.Kn[3], 0, 0.01, 0.0, 0.001 };
	variables[ 8] = (optimum_variable){ AdjustGlobal, &camcal.Kn[4], 0, 0.01, 0.0, 0.001 };
	variables[ 9] = (optimum_variable){ AdjustGlobal, &camcal.Kn[5], 0, 0.01, 0.0, 0.001 };
	variables[10] = (optimum_variable){ AdjustGlobal, &camcal.Pn[0], 0, 0.01, 0.0, 0.001 };
	variables[11] = (optimum_variable){ AdjustGlobal, &camcal.Pn[1], 0, 0.01, 0.0, 0.001 };
	variables[12] = (optimum_variable){ AdjustGlobal, &camcal.Sn[0], 0, 0.01, 0.0, 0.001 };
	variables[13] = (optimum_variable){ AdjustGlobal, &camcal.Sn[1], 0, 0.01, 0.0, 0.001 };
	variables[14] = (optimum_variable){ AdjustGlobal, &camcal.Sn[2], 0, 0.01, 0.0, 0.001 };
	variables[15] = (optimum_variable){ AdjustGlobal, &camcal.Sn[3], 0, 0.01, 0.0, 0.001 };

	optimum_ctx ctx = {
		.fnCompute = GetGlobalError,
		.fnContinueGoing = ContinueGoing,
		.vOpaque = 0,
		.iOpaque = 0,
		.nrVariables = tuningVars,
		.variables = variables,
	};
	double er = SFHGenericOptimumOptimize( &ctx );
	return er;
}

int main()
{
	SetCalibrationDefaults( &camcal, calimages, nrcalimages );

	double err = Calibrate( &camcal, calimages, nrcalimages );

	printf( "Final error: %f\n", err );
	if( err > 0.002 )
	{
		fprintf( stderr, "FAIL: Error too high\n" );
		return -1;
	}
	else
	{
		fprintf( stderr, "PASS: Error within tolerance\n" );
		return 0;
	}
}


