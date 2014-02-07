#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include "wells.h"
#include "design.h"
// LAPLACE INVERSION VARIABLES
//double DHTOL = 1.0e-6,DHALPHA=0.0,time_max = 10.2, tol = 1.0e-6;

double inversion( double _Complex *conc_lap, double time, double max_time, double DHTOL, double DHALPHA, int max_lap );
extern void __complex_bessel_MOD_cbesj( double _Complex *, double *, int *, int *, double _Complex *, int *, int * );
extern void __complex_bessel_MOD_cbesk( double _Complex *, double *, int *, int *, double _Complex *, int *, int * );
extern void __complex_bessel_MOD_cbesy( double _Complex *, double *, int *, int *, double _Complex *, int *, int * );
double complex unc_lap_sat( double complex, double y, double ts, struct params *p );
double complex leaky_unc_lap_sat( double complex, double y, double ts, struct params *p );
double complex sC_zD( double complex, double y, double ts, double zD, struct params *p );
double complex confined_lap( double complex p_i, double ts, struct params *p );
double BESSJ0( double X );
double BESSJ1( double X );

double unc_Step_Lap_sat( double y, void *params )
{
//double PI = 3.142857142857142793701541449991054832935333251953125;
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 10.2,p0,tee;
	struct params *p = (struct params * )params;
	double  time_inv = 1.0, p0, tee;
	int i, max_lap = 40;
	double complex lap_func[max_lap], p_i;
	double sU;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	for( i = 0; i < max_lap; i++ )
	{
		p_i = p0 + ( double )( i ) * ( M_PI / tee ) * I;
		lap_func[i] = unc_lap_sat( p_i, y, ( *p ).ts, p );
	}
	sU = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap ) * y * BESSJ0( sqrt( ( *p ).beta2 ) * y );
	return( sU );
}

double unc_Linear_Lap_sat( double y, void *params )
{
//double PI = 3.142857142857142793701541449991054832935333251953125;
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 10.2,p0,tee;
	struct params *p = (struct params * )params;
	double  time_inv = 1.0, p0, tee;
	int i, nLoop, j, max_lap = 40;
	double complex lap_func[max_lap], p_i;
	double sd[2], td[2], Q[2], beta;
	double sU;
	beta = ( ( *p ).Q2 - ( *p ).Q1 ) / ( ( *p ).ts2 - ( *p ).ts1 );
	Q[0] = ( *p ).Q1; Q[1] = ( *p ).Q2;
	td[0] = ( *p ).ts - ( *p ).ts1; td[1] = ( *p ).ts - ( *p ).ts2;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	if( ( *p ).ts <= ( *p ).ts1 ) {sU = 0.0; nLoop = 0; return sU; exit( 0 );}
	if( ( *p ).ts <= ( *p ).ts2 ) nLoop = 1;
	else nLoop = 2;
	sd[0] = sd[1] = 0.0;
	for( j = 0; j < nLoop; j++ ) // j = 0, 1
	{
		for( i = 0; i < max_lap; i++ )
		{
			p_i = p0 + ( double )( M_PI * i / tee ) * I;
			lap_func[i] = ( double complex )( 1 - 2 * j ) * ( Q[j] + beta * td[j] / p_i ) * unc_lap_sat( p_i, y, td[j], p );
		}
		sd[j] = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
	}
	sU = ( sd[0] + sd[1] ) * y * BESSJ0( sqrt( ( *p ).beta2 ) * y );
	return( sU );
}

double complex unc_lap_sat( double complex p_i, double y, double ts, struct params *p )
{
	double complex mu, musq, qD_inv, lap_sat;
	double complex lambda, lambda_1, lambda_2, chai, B_D;
	double complex phi_zero, phi_L, ctemp1, ctemp, bj[4], by[4];
	int nz = 0, ierr = 0, n = 1, kode = 1;
	double fnuN, fnuM;
	double akD = ( *p ).akD;
	double acD = ( *p ).acD;
	double l_unsat = ( *p ).l_unsat;
	musq = y * y + p_i / ( ts * ( *p ).beta2 );
	mu = csqrt( musq );
	B_D = p_i * acD * exp( akD * ( *p ).psiD ) / ( ( *p ).sigma * ts * ( *p ).beta2 );
	if( abs( akD - acD ) == 0.0 )
	{
		lambda = csqrt( B_D + y * y + akD * akD / ( double )4.0 );
		lambda_1 = -lambda + akD / ( double )2.0;
		lambda_2 = lambda + akD / ( double )2.0;
		chai = -( lambda_1 / lambda_2 ) * cexp( -lambda * l_unsat * 2.0 );
		qD_inv = ( double complex )( chai + 1.0 ) / ( lambda_1 + chai * lambda_2 );
//printf("qD_inv = %g %g\n",qD_inv);
	}
	else
	{
		fnuN = sqrt( ( akD * akD + y * y * 4.0 ) / ( ( akD - acD ) * ( akD - acD ) ) );
		fnuM = fnuN + 1.0;
// Compute phi_0 at water Table
		phi_zero = ( double complex )( 0.0 + 1.0 * I ) * csqrt( B_D * 4.0 / ( ( akD - acD ) * ( akD - acD ) ) );
		__complex_bessel_MOD_cbesj( &phi_zero, &fnuN, &n, &kode, &bj[0], &nz, &ierr ) ; //Jn(phi_zero)
		__complex_bessel_MOD_cbesy( &phi_zero, &fnuN, &n, &kode, &by[0], &nz, &ierr ) ; //Yn(phi_zero)
		__complex_bessel_MOD_cbesj( &phi_zero, &fnuM, &n, &kode, &bj[1], &nz, &ierr ) ; //Jm(phi_zero)
		__complex_bessel_MOD_cbesy( &phi_zero, &fnuM, &n, &kode, &by[1], &nz, &ierr ) ; //Ym(phi_zero)
// Compute chai for small or Large Arguments
		if( ( akD * l_unsat ) > 1.0 )
		{
			chai = 0.0 + 1.0 * I;
		}
		else
		{
			phi_L = ( 0.0 + 1.0 * I ) * csqrt( B_D * ( double )4.0 * exp( ( akD - acD ) * l_unsat ) / ( ( akD - acD ) * ( akD - acD ) ) );
			__complex_bessel_MOD_cbesj( &phi_L, &fnuN, &n, &kode, &bj[2], &nz, &ierr ) ; //Jn(phi_zero)
			__complex_bessel_MOD_cbesy( &phi_L, &fnuN, &n, &kode, &by[2], &nz, &ierr ) ; //Yn(phi_zero)
			__complex_bessel_MOD_cbesj( &phi_L, &fnuM, &n, &kode, &bj[3], &nz, &ierr ) ; //Jm(phi_zero)
			__complex_bessel_MOD_cbesy( &phi_L, &fnuM, &n, &kode, &by[3], &nz, &ierr ) ; //Ym(phi_zero)
			chai = akD * ( double )( 0.5 ) + sqrt( akD * akD * ( double )( 0.25 ) + y * y ) * bj[2] - ( 0.0 + 1.0 * I ) * csqrt( B_D * exp( akD * l_unsat ) ) * bj[3];
			chai = chai / akD * ( double )( 0.5 ) + sqrt( akD * akD * ( double )( 0.25 ) + y * y ) * by[2] - ( 0.0 + 1.0 * I ) * csqrt( B_D * exp( akD * l_unsat ) ) * by[3];
		}
// Compute qD with caution
// Check denominator of qD does not goes to zero
		ctemp1 = bj[1] + chai * by[1];
		ctemp  = bj[0] + chai * by[0];
		if( abs( ctemp ) == abs( ctemp1 ) )
		{
			qD_inv = ( double complex ) 1.0 / ( akD * 0.5 + sqrt( 0.25 * akD * akD + y * y ) - ( 0.0 + 1.0 * I ) * csqrt( B_D ) );
		}
		else if( abs( ctemp ) == 0.0 )
		{
			qD_inv = 0.0 + 0.0 * I;
		}
		else
		{
			qD_inv = ( double complex ) 1.0 / ( akD * 0.5 + sqrt( 0.25 * akD * akD + y * y ) - ( 0.0 + 1.0 * I ) * csqrt( B_D ) * ctemp1 / ctemp );
		}
	}
// Now Compute the Final form of Laplace Transformed sU
	if( ( *p ).zD1 == ( *p ).zD2 )
	{
		lap_sat = sC_zD( p_i, y, ts, 1.0, p ) * ccosh( mu * ( ( double )1.0 - ( *p ).zD1 ) ) / ( ccosh( mu ) - ( mu * qD_inv ) * csinh( mu ) );
	}
	else
	{
		lap_sat = ( sC_zD( p_i, y, ts, 1.0, p ) / ( mu * ( ( *p ).zD2 - ( *p ).zD1 ) ) ) * ( csinh( mu * ( ( double )1.0 - ( *p ).zD1 ) ) - csinh( mu * ( ( double )1.0 - ( *p ).zD2 ) ) ) / ( ccosh( mu ) - ( mu * qD_inv ) * csinh( mu ) );
	}
//printf("akD = %g \n", (*p).acD);
	return ( lap_sat );
}

/***********************************/
/*                                 */
/*  LEAKY UNCONFINED AQUIFER PART  */
/*                                 */
/***********************************/

double leaky_unc_Step_Lap_sat( double y, void *params )
{
//double PI = 3.142857142857142793701541449991054832935333251953125;
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 10.2,p0,tee;
	struct params *p = (struct params * )params;
	double  time_inv = 1.0, p0, tee;
	int i, max_lap = 40;
	double complex lap_func[max_lap], p_i;
	double sU;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	for( i = 0; i < max_lap; i++ )
	{
		p_i = p0 + ( double )( i ) * ( M_PI / tee ) * I;
		lap_func[i] = leaky_unc_lap_sat( p_i, y, ( *p ).ts, p );
	}
	sU = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap ) * y * BESSJ0( sqrt( ( *p ).beta2 ) * y );
	return( sU );
}

double leaky_unc_Linear_Lap_sat( double y, void *params )
{
//double PI = 3.142857142857142793701541449991054832935333251953125;
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 10.2,p0,tee;
	struct params *p = (struct params * )params;
	double  time_inv = 1.0, p0, tee;
	int i, nLoop, j, max_lap = 40;
	double complex lap_func[max_lap], p_i;
	double sd[2], td[2], Q[2], beta;
	double sU;
	beta = ( ( *p ).Q2 - ( *p ).Q1 ) / ( ( *p ).ts2 - ( *p ).ts1 );
	Q[0] = ( *p ).Q1; Q[1] = ( *p ).Q2;
	td[0] = ( *p ).ts - ( *p ).ts1; td[1] = ( *p ).ts - ( *p ).ts2;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	if( ( *p ).ts <= ( *p ).ts1 ) {sU = 0.0; nLoop = 0; return sU; exit( 0 );}
	if( ( *p ).ts <= ( *p ).ts2 ) nLoop = 1;
	else nLoop = 2;
	sd[0] = sd[1] = 0.0;
	for( j = 0; j < nLoop; j++ ) // j = 0, 1
	{
		for( i = 0; i < max_lap; i++ )
		{
			p_i = p0 + ( double )( M_PI * i / tee ) * I;
			lap_func[i] = ( double complex )( 1 - 2 * j ) * ( Q[j] + beta * td[j] / p_i ) * leaky_unc_lap_sat( p_i, y, td[j], p );
		}
		sd[j] = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
	}
	sU = ( sd[0] + sd[1] ) * y * BESSJ0( sqrt( ( *p ).beta2 ) * y );
	return( sU );
}

double complex leaky_unc_lap_sat( double complex p_i, double y, double ts, struct params *p )
{
	double complex mu, musq, qD_inv, lap_sat;
	double complex lambda, lambda_1, lambda_2, chai, B_D;
	double complex phi_zero, phi_L, ctemp1, ctemp, bj[4], by[4];
	int nz = 0, ierr = 0, n = 1, kode = 1;
	double fnuN, fnuM;
	double akD = ( *p ).akD;
	double acD = ( *p ).acD;
	double l_unsat = ( *p ).l_unsat;
	musq = y * y + p_i / ( ts * ( *p ).beta2 );
	mu = csqrt( musq );
	B_D = p_i * acD * exp( akD * ( *p ).psiD ) / ( ( *p ).sigma * ts * ( *p ).beta2 );
	if( abs( akD - acD ) == 0.0 )
	{
		lambda = csqrt( B_D + y * y + akD * akD / ( double )4.0 );
		lambda_1 = -lambda + akD / ( double )2.0;
		lambda_2 = lambda + akD / ( double )2.0;
		chai = -( lambda_1 / lambda_2 ) * cexp( -lambda * l_unsat * 2.0 );
		qD_inv = ( double complex )( chai + 1.0 ) / ( lambda_1 + chai * lambda_2 );
//printf("qD_inv = %g %g\n",qD_inv);
	}
	else
	{
		fnuN = sqrt( ( akD * akD + y * y * 4.0 ) / ( ( akD - acD ) * ( akD - acD ) ) );
		fnuM = fnuN + 1.0;
// Compute phi_0 at water Table
		phi_zero = ( double complex )( 0.0 + 1.0 * I ) * csqrt( B_D * 4.0 / ( ( akD - acD ) * ( akD - acD ) ) );
		__complex_bessel_MOD_cbesj( &phi_zero, &fnuN, &n, &kode, &bj[0], &nz, &ierr ) ; //Jn(phi_zero)
		__complex_bessel_MOD_cbesy( &phi_zero, &fnuN, &n, &kode, &by[0], &nz, &ierr ) ; //Yn(phi_zero)
		__complex_bessel_MOD_cbesj( &phi_zero, &fnuM, &n, &kode, &bj[1], &nz, &ierr ) ; //Jm(phi_zero)
		__complex_bessel_MOD_cbesy( &phi_zero, &fnuM, &n, &kode, &by[1], &nz, &ierr ) ; //Ym(phi_zero)
// Compute chai for small or Large Arguments
		if( ( akD * l_unsat ) > 1.0 )
		{
			chai = 0.0 + 1.0 * I;
		}
		else
		{
			phi_L = ( 0.0 + 1.0 * I ) * csqrt( B_D * ( double )4.0 * exp( ( akD - acD ) * l_unsat ) / ( ( akD - acD ) * ( akD - acD ) ) );
			__complex_bessel_MOD_cbesj( &phi_L, &fnuN, &n, &kode, &bj[2], &nz, &ierr ) ; //Jn(phi_zero)
			__complex_bessel_MOD_cbesy( &phi_L, &fnuN, &n, &kode, &by[2], &nz, &ierr ) ; //Yn(phi_zero)
			__complex_bessel_MOD_cbesj( &phi_L, &fnuM, &n, &kode, &bj[3], &nz, &ierr ) ; //Jm(phi_zero)
			__complex_bessel_MOD_cbesy( &phi_L, &fnuM, &n, &kode, &by[3], &nz, &ierr ) ; //Ym(phi_zero)
			chai = akD * ( double )( 0.5 ) + sqrt( akD * akD * ( double )( 0.25 ) + y * y ) * bj[2] - ( 0.0 + 1.0 * I ) * csqrt( B_D * exp( akD * l_unsat ) ) * bj[3];
			chai = chai / akD * ( double )( 0.5 ) + sqrt( akD * akD * ( double )( 0.25 ) + y * y ) * by[2] - ( 0.0 + 1.0 * I ) * csqrt( B_D * exp( akD * l_unsat ) ) * by[3];
		}
// Compute qD with caution
// Check denominator of qD does not goes to zero
		ctemp1 = bj[1] + chai * by[1];
		ctemp  = bj[0] + chai * by[0];
		if( abs( ctemp ) == abs( ctemp1 ) )
		{
			qD_inv = ( double complex ) 1.0 / ( akD * 0.5 + sqrt( 0.25 * akD * akD + y * y ) - ( 0.0 + 1.0 * I ) * csqrt( B_D ) );
		}
		else if( abs( ctemp ) == 0.0 )
		{
			qD_inv = 0.0 + 0.0 * I;
		}
		else
		{
			qD_inv = ( double complex ) 1.0 / ( akD * 0.5 + sqrt( 0.25 * akD * akD + y * y ) - ( 0.0 + 1.0 * I ) * csqrt( B_D ) * ctemp1 / ctemp );
		}
	}
// Now Compute the Final form of Laplace Transformed sU
	if( ( *p ).zD1 == ( *p ).zD2 )
	{
		lap_sat = sC_zD( p_i, y, ts, 1.0, p ) * ccosh( mu * ( ( double )1.0 - ( *p ).zD1 ) ) / ( ccosh( mu ) - ( mu * qD_inv ) * csinh( mu ) );
	}
	else
	{
		lap_sat = ( sC_zD( p_i, y, ts, 1.0, p ) / ( mu * ( ( *p ).zD2 - ( *p ).zD1 ) ) ) * ( csinh( mu * ( ( double )1.0 - ( *p ).zD1 ) ) - csinh( mu * ( ( double )1.0 - ( *p ).zD2 ) ) ) / ( ccosh( mu ) - ( mu * qD_inv ) * csinh( mu ) );
	}
//printf("akD = %g \n", (*p).acD);
	return ( lap_sat );
}

double complex sC_zD( double complex p_i, double y, double ts, double zD, struct params *p )
{
	double complex sC, musq, mu, z1, z2, bj[2], bk[4], phi_0, phi_n, ctemp, ctemp1, fac2, fac1;
	double PI2, yD, beta, fnu0 = 0.0, fnu1 = 1.0;
	int nz = 0, ierr = 0, k = 1, kode = 1;
	int n;
	double dD = ( *p ).dD;
	double lD = ( *p ).lD;
	double rwD = ( *p ).rwD;
	double cwD = ( *p ).cwD;
	double beta2 = ( *p ).beta2;
	musq = y * y + p_i / ( ts * beta2 );
	mu = csqrt( musq );
//printf("rwD=%g \n",rwD);
	if( rwD == 0.0 )
	{
		sC = ( double complex )( -2.0 / p_i ) * ( csinh( mu * ( ( double )1.0 - lD ) ) - csinh( mu * ( ( double )1.0 - dD ) ) ) / ( musq * ( lD - dD ) * csinh( mu ) );
	}
	else
	{
		PI2 = M_PI * M_PI;
		beta = sqrt( beta2 );
		yD = y * beta;
		phi_0 = csqrt( p_i / ts );
		z1 = rwD * phi_0;
		__complex_bessel_MOD_cbesk( &z1, &fnu0, &k, &kode, &bk[0], &nz, &ierr ) ; //K0(phi_)*rwD)
		__complex_bessel_MOD_cbesk( &z1, &fnu1, &k, &kode, &bk[1], &nz, &ierr ) ; //K0(phi_)*rwD)
		bj[0] = BESSJ0( yD * rwD ); bj[1] = BESSJ1( yD * rwD );
		n = 1;
		ctemp = 0.0 + 0.0 * I;
		do
		{
			phi_n = csqrt( p_i / ts + beta2 * PI2 * n * n );
			z2 = rwD * phi_n;
			__complex_bessel_MOD_cbesk( &z2, &fnu0, &k, &kode, &bk[2], &nz, &ierr ) ; //K0(phi_n*rwD)
			__complex_bessel_MOD_cbesk( &z2, &fnu1, &k, &kode, &bk[3], &nz, &ierr ) ; //K0(phi_n*rwD)
			ctemp1 = ( double complex )( ( sin( n * M_PI * lD ) - sin( n * M_PI * dD ) ) / ( phi_n * bk[3] + 0.5 * cwD * rwD * phi_0 * phi_0 * bk[2] / ( lD - dD ) ) );
			if( y == 0.0 )
			{
				fac2 = phi_n * bk[3] / ( musq + n * n * PI2 );
			}
			else
			{
				fac2 = ( double complex )( ( beta / y ) * bj[1] * bk[2] + ( phi_n * bj[0] * bk[3] - yD * bj[1] * bk[2] ) / ( musq + n * n * PI2 ) );
			}
			//ctemp1=ctemp1*fac2/n;
			ctemp1 = ctemp1 * fac2 * cos( n * M_PI * ( 1.0 - zD ) ) / n;
			//if(isnan(abs(ctemp1))) break;
			ctemp = ctemp + ctemp1;
			n = n + 1;
		}
		while( cabs( ctemp1 ) > ( *p ).DHTOL );
		if( y == 0.0 )
		{
			fac1 = phi_0 * bk[1] / musq;
		}
		else
		{
			fac1 = ( double complex )( ( beta / y ) * bj[1] * bk[0] + ( phi_0 * bj[0] * bk[1] - yD * bj[1] * bk[0] ) / musq );
		}
		sC = ( double complex )( ctemp * 4.0 / ( p_i * M_PI * ( lD - dD ) ) + ( 2.0 / p_i ) * fac1 / ( phi_0 * bk[1] + 0.5 * cwD * rwD * phi_0 * phi_0 * bk[0] / ( lD - dD ) ) );
	}
	sC = -sC;
//printf("sC_zD=%g \n",sC);
	return ( sC );
}


double confined_Step_Lap( struct params *p )
{
	/* This Function computes dimesnionless Drawdwon at Dimensionless Time ts */
	/* ts=(Kr/Ss)*t/r2 dim_draw=4*pi*T*s/Q */
	double  time_inv = 1.0, p0, tee, dim_draw;
	int max_lap = 40, i;
	double complex lap_func[max_lap], p_i;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	for( i = 0; i < max_lap; i++ )
	{
		p_i = p0 + ( ( double )( i ) * ( M_PI / tee ) ) * I;
		lap_func[i] = ( double complex ) confined_lap( p_i, ( *p ).ts, p );
	}
	dim_draw = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
	return dim_draw;
}

/*User Defined lap ***/
double confined_user_def_Lap( struct params *p, struct Well_Data *wd, int k , double tf )
{
	double  time_inv = 1.0, p0, tee;
	int max_lap = 40, j, nR;
	double complex lap_func[max_lap], p_i, Qp;
	double dim_draw;
	nR = ( *wd ).nQ[k];
	dim_draw = 0.0;
//printf ("nR = %d \n",nR);
//tf = tf/(*p).ts;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	Qp = 0.0 + 0.0 * I;
	for( j = 0; j < max_lap; j++ )
	{
		p_i = p0 + ( double )( M_PI * j / tee ) * I;
		/*step-shutoff */
		// Qp = (2.0)*(1.0-cexp(-0.5*tf*p_i/(*p).ts));
		/* Sinusoidal with sin(wt) +C form */
		Qp = 2.0 + ( ( 30 / M_PI ) ) / ( ( 30 / M_PI ) * ( 30 / M_PI ) * ( *p ).ts / ( tf * p_i ) + p_i * tf / ( *p ).ts );
		lap_func[j] = ( Qp ) * confined_lap( p_i, ( *p ).ts, p );
		//printf("Lap = %g %g  \n",j, Qp2);
	}
	dim_draw = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
//  printf("dim_draw = %g \n", dim_draw);
	return dim_draw;
}
/*Cumulative step ***/
double confined_Step_cumulative_Lap( struct params *p, struct Well_Data *wd, int k , double tf )
{
	double  time_inv = 1.0, p0, tee;
	int max_lap = 40, i, j, nR;
	double complex lap_func[max_lap], p_i, Qp;
	double dim_draw;
	nR = ( *wd ).nQ[k];
	dim_draw = 0.0;
//printf ("nR = %d \n",nR);
//tf = tf/(*p).ts;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	Qp = 0.0 + 0.0 * I;
	for( j = 0; j < max_lap; j++ )
	{
		p_i = p0 + ( double )( M_PI * j / tee ) * I;
		Qp = ( *wd ).Q[k][0] + 0.0 * I;
		for( i = 1; i < nR && ( ( *wd ).t[k][i - 1] - ( *wd ).t[k][0] )*tf  < ( *p ).ts; i++ )
		{
			Qp += cexp( -p_i * ( ( *wd ).t[k][i] - ( *wd ).t[k][0] ) * tf / ( *p ).ts ) * ( ( *wd ).Q[k][i] - ( *wd ).Q[k][i - 1] );
		}
		lap_func[j] = ( Qp ) * confined_lap( p_i, ( *p ).ts, p );
		//printf("Lap = %g %g  \n",j, Qp2);
	}
	dim_draw = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
//  printf("dim_draw = %g \n", dim_draw);
	return dim_draw;
}
//***************
double confined_Linear_convolute_Lap( struct params *p )
{
	/* This Function computes dimesnionless Drawdwon at Dimensionless Time ts */
// USING CONVOLUTION OF PUMPING DUE TO LINEAR ELEMENTS
	/* ts=(Kr/Ss)*t/r2 dim_draw=4*pi*T*s/Q */
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 1.2,p0,tee;
	double  time_inv = 1.0, p0, tee;
	int max_lap = 40, i, j, nLoop;
	double complex lap_func[max_lap], p_i;
	double sd[2], dim_draw, td[2], Q[2], beta;
	beta = ( ( *p ).Q2 - ( *p ).Q1 ) / ( ( *p ).ts2 - ( *p ).ts1 );
	Q[0] = ( *p ).Q1; Q[1] = ( *p ).Q2;
	td[0] = ( *p ).ts - ( *p ).ts1; td[1] = ( *p ).ts - ( *p ).ts2;
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	if( ( *p ).ts <= ( *p ).ts1 ) {dim_draw = 0.0; nLoop = 0; return dim_draw; exit( 0 );}
	if( ( *p ).ts <= ( *p ).ts2 ) nLoop = 1;
	else nLoop = 2;
	sd[0] = sd[1] = 0.0;
	for( j = 0; j < nLoop; j++ ) // j = 0, 1
	{
		for( i = 0; i < max_lap; i++ )
		{
			p_i = p0 + ( double )( M_PI * i / tee ) * I;
			lap_func[i] = ( double complex )( 1 - 2 * j ) * ( Q[j] + beta * td[j] / p_i ) * confined_lap( p_i, td[j], p );
		}
		sd[j] = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
	}
	dim_draw = sd[0] + sd[1];
	return dim_draw;
}

double confined_Linear_Lap( struct params *p, struct Well_Data *wd, int k, double tf )
{
	/* This Function computes dimesnionless Drawdwon at Dimensionless Time ts */
// USING DIRECT METHOD
	/* ts=(Kr/Ss)*t/r2 dim_draw=4*pi*T*s/Q */
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 1.2,p0,tee;
	double  time_inv = 1.0, p0, tee;
	int max_lap = 40, i, j, nR;
	double complex lap_func[max_lap], p_i, Qp1, Qp2;
	double dim_draw, beta;
	nR = ( *wd ).nQ[k];
//printf ("nR = %d \n",nR);
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	for( j = 0; j < max_lap; j++ )
	{
		p_i = p0 + ( double )( M_PI * j / tee ) * I;
		Qp1 = 0.0 + 0.0 * I; Qp2 = 0.0 + 0.0 * I;
		for( i = 1; i < nR ; i++ )
		{
			beta = ( ( *wd ).Q[k][i] - ( *wd ).Q[k][i - 1] ) / ( ( *wd ).t[k][i] * tf - ( *wd ).t[k][i - 1] * tf );
			Qp1 += ( ( *wd ).Q[k][i - 1] + beta * ( *p ).ts / p_i ) * ( cexp( -( *wd ).t[k][i - 1] * p_i * tf / ( *p ).ts ) - cexp( -( *wd ).t[k][i] * p_i * tf / ( *p ).ts ) );
			Qp2 += beta * ( ( *wd ).t[k][i] * tf - ( *wd ).t[k][i - 1] * tf ) * cexp( -( *wd ).t[k][i] * p_i * tf / ( *p ).ts );
//       printf("Qp1 = %g %g Qp2 = %g %g \n", Qp1, Qp2);
//       printf("beta = %g  \n", beta);
		}
		lap_func[j] = ( Qp1 - Qp2 ) * confined_lap( p_i, ( *p ).ts, p );
		//printf("Lap = %g %g  \n",j, Qp2);
	}
	dim_draw = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
//  printf("dim_draw = %g \n", dim_draw);
	return dim_draw;
}

double confined_Linear_superpose_Lap( struct params *p, struct Well_Data *wd, int k, double tf )
{
	/* This Function computes dimesnionless Drawdwon at Dimensionless Time ts */
// USING SUPERPOSITION OF LINEAR ELEMENTS
	/* ts=(Kr/Ss)*t/r2 dim_draw=4*pi*T*s/Q */
//double  DHTOL =1.0e-6,DHALPHA=0.0,time_inv = 1.0,time_max = 1.2,p0,tee;
	double  time_inv = 1.0, p0, tee;
	int max_lap = 40, i, j, nR;
	double complex lap_func[max_lap], p_i, Qp1, Qp2;
	double dim_draw, beta;
	nR = ( *wd ).nQ[k];
	tee = ( double ) 0.8 * ( *p ).time_max;
	p0 = ( ( *p ).DHALPHA - log( ( *p ).DHTOL ) ) / ( ( double ) 2.00 * tee );
	for( j = 0; j < max_lap; j++ )
	{
		p_i = p0 + ( double )( M_PI * j / tee ) * I;
		Qp1 = 0.0 + 0.0 * I; Qp2 = 0.0 + 0.0 * I;
		for( i = 1; i < nR ; i++ )
		{
			beta = ( ( *wd ).Q[k][i] - ( *wd ).Q[k][i - 1] ) / ( ( *wd ).t[k][i] * tf - ( *wd ).t[k][i - 1] * tf );
			Qp1 += ( ( *wd ).Q[k][i - 1] + beta * ( *p ).ts / p_i ) * ( cexp( -( *wd ).t[k][i - 1] * p_i * tf / ( *p ).ts ) );
			Qp2 += ( ( *wd ).Q[k][i] + beta * ( *p ).ts / p_i ) * ( cexp( -( *wd ).t[k][i] * p_i * tf / ( *p ).ts ) );
		}
		lap_func[j] = ( Qp1 - Qp2 ) * confined_lap( p_i, ( *p ).ts, p );
	}
	dim_draw = inversion( lap_func, time_inv, ( *p ).time_max, ( *p ).DHTOL, ( *p ).DHALPHA, max_lap );
	if( isnan( dim_draw ) ) dim_draw = 1.0e-10;
	return dim_draw;
}

/* Function To get Laplace Transformed Drawdown in Confined Aquifer Mishra and Neuman (2011) */
double complex confined_lap( double complex p_i, double ts, struct params *p )
{
	/* This Function computes dimesnionless Drawdwon at Dimensionless Time ts */
	/* ts=(Kr/Ss)*t/r2 dim_draw=4*pi*T*s/Q */
	double complex lap_func, ctemp, ctemp1, s0;
	double complex z0, z1, z2, z3, denom1, denom2, bk0[4], bk1[2];
	int nz = 0, ierr = 0, n = 1, kode = 1, k;
	double fnu0 = 0.0, fnu1 = 1.0, s1;
	double dD = ( *p ).dD;
	double lD = ( *p ).lD;
	double rwD = ( *p ).rwD;
//printf("%g %g \n", (*p).zD1,(*p).zD2);
	z0 = csqrt( p_i / ts );
	__complex_bessel_MOD_cbesk( &z0, &fnu0, &n, &kode, &bk0[0], &nz, &ierr );
	denom1 = 1.0; denom2 = 1.0;
	if( rwD > 0 )
	{
		z1 = rwD * z0;
		__complex_bessel_MOD_cbesk( &z1, &fnu0, &n, &kode, &bk0[1], &nz, &ierr );
		__complex_bessel_MOD_cbesk( &z1, &fnu1, &n, &kode, &bk1[0], &nz, &ierr );
		denom1 = ( double complex )( z1 * bk1[0] + 0.5 * ( *p ).cwD * z1 * z1 * bk0[1] / ( lD - dD ) );
	}
	s0 = p_i / ts;
	s1 = ( *p ).beta2 * M_PI * M_PI;
	ctemp = 0.0 + 0.0 * I; ctemp1 = 0.0 + 0.0 * I; k = 0;
	if( ( *p ).lD != 1.0 && ( *p ).dD != 0.0 )
	{
		do
		{
			k++;
			z2 = csqrt( s0 + s1 * k * k );
			__complex_bessel_MOD_cbesk( &z2, &fnu0, &n, &kode, &bk0[2], &nz, &ierr ) ;
			if( rwD > 0 )
			{
				z3 = rwD * z2;
				__complex_bessel_MOD_cbesk( &z3, &fnu0, &n, &kode, &bk0[3], &nz, &ierr ) ;
				__complex_bessel_MOD_cbesk( &z3, &fnu1, &n, &kode, &bk1[1], &nz, &ierr ) ;
				denom2 = ( double complex )( z3 * bk1[1] + 0.5 * ( *p ).cwD * z1 * z1 * +bk0[3] / ( lD - dD ) );
			}
			/* Check if point or average drawdown needed */
			/* zD1 and zD2 are depths below the top boundary */
			if( ( *p ).zD1 == ( *p ).zD2 )
			{
				ctemp += ctemp1 = ( double complex ) bk0[2] * ( sin( M_PI * lD * k ) - sin( M_PI * dD * k ) ) * cos( M_PI * ( *p ).zD1 * k ) / ( k * denom2 ) ;
			}
			else
			{
				ctemp += ctemp1 = ( double complex )( bk0[2] / ( M_PI * ( ( *p ).zD2 - ( *p ).zD1 ) ) ) * ( sin( M_PI * lD * k ) - sin( M_PI * dD * k ) ) * ( sin( M_PI * ( *p ).zD2 * k ) - sin( M_PI * ( *p ).zD1 * k ) ) / ( k * denom2 ) ;
			}
		}
		while( cabs( ctemp1 ) > ( *p ).DHTOL );
	}
	lap_func = ( double complex ) 2.0 * bk0[0] / ( denom1 * p_i ) + ctemp * 4.0 / ( p_i * M_PI * ( lD - dD ) );
	return lap_func;
}


/***************************-------------------**********************************************/
/**************************| LAPLACE INVERSION |********************************************/
/****************************------------------*********************************************/
void  add( int count, double complex *conc, double _Complex **epsilon, int *traffic );
void am_bm( double complex *record, double complex z, double complex *am, double complex *bm, int max_lap );
void coeff_QD( double complex *record, double complex *conc, int max_lap );
double complex Term( int l, double complex *conc, int *traffic );
double inversion( double complex *conc_lap, double time, double max_time, double DHTOL, double DHALPHA, int max_lap )
{
	int i;
	double p0, tee, ret_value;
	double complex record[max_lap], conc[max_lap], temp_time, am, bm;
	tee = ( double ) 0.8 * max_time;
	p0 = ( DHALPHA - log( DHTOL ) ) / ( ( double ) 2.00 * tee );
	/* Initialize conc and record */
	for( i = 0; i < max_lap; i++ )
	{
		conc[i] = conc_lap[i];
		record[i] = ( double )0.0 + ( double )0.0 * I;
	}
	temp_time = cexp( ( double )0.0 + ( M_PI * time / ( ( double )0.8 * max_time ) ) * I );
	am = 0.0 + 0.0 * I; bm = 0.0 + 0.0 * I;
	coeff_QD( record, conc, max_lap );
	am_bm( record, temp_time, &am, &bm, max_lap );
	if( cabs( am ) == 0.0 )
	{ ret_value = ( double )0.0;}
	else if( cabs( bm ) == 0.0 )
	{ printf( "%s", "LAPLACE INVERSION ERROR" ); exit( -1 );}
	else
	{ ret_value = cexp( p0 * time ) / ( double )0.8 / max_time * creal( am / bm );}
	return ret_value;
}

/*Function Coeff_QD */

void coeff_QD( double complex *record, double complex *conc, int max_lap )
{
	int i, j, count = 0, traffic = 1;
	double complex **epsilon;
	epsilon = ( double complex ** ) malloc( max_lap * sizeof( double complex * ) );
	for( i = 0; i < max_lap; i++ )
		epsilon[i] = ( double complex * ) malloc( 2 * sizeof( double complex ) );
	/* Initialize epsilon array*/
	for( i = 0; i < 2; i++ )
	{
		for( j = 0; j < max_lap; j++ )
		{epsilon[j][i] = ( double )0.0 + ( double )0.0 * I;}
	}
	record[count] = Term( -1, conc, &traffic );
	for( i = 0; i < max_lap - 1; i++ )
	{
		add( i, conc, epsilon, &traffic );
		count = count + 1;
		record[count] = -epsilon[0][0];
	}
}

/*Function ab to get am bm */
void am_bm( double complex *record, double complex z, double complex *am, double complex *bm, int max_lap )
{
	int i;
	double complex temp_1[max_lap], temp_2[max_lap], h;
	temp_1[0] = record[0];
	temp_1[1] = record[0];
	temp_2[0] = ( double )1.0 + ( double )0.0 * I;
	temp_2[1] = ( double )1.0 + record[1] * z;
	for( i = 2; i < max_lap; i++ )
	{
		if( i < ( max_lap - 1 ) )
		{
			temp_1[i] = temp_1[i - 1] + record[i] * temp_1[i - 2] * z;
			temp_2[i] = temp_2[i - 1] + record[i] * temp_2[i - 2] * z;
		}
		else if( i == ( max_lap - 1 ) )
		{
			h = ( double )0.5 + ( double )0.5 * z * ( record[i - 1] - record[i] );
			h = -h * ( ( double )1 - cpow( ( ( double )1 + z * record[i] / ( h * h ) ), ( double )0.5 ) );
			temp_1[i] = temp_1[i - 1] + h * temp_1[i - 2] * z;
			temp_2[i] = temp_2[i - 1] + h * temp_2[i - 2] * z;
			*am = temp_1[i];
			*bm = temp_2[i];
		}
	}
}

/*Function Term */
double complex Term( int l, double complex *conc, int *traffic )
{
	double complex ret_value ;
	ret_value = conc[l + 1];
	*traffic = 1;
	if( l < 0 ) ret_value = ( ( double )0.5 ) * conc[0];
	if( cabs( ret_value ) < exp( ( double ) - 1.0 * 50.0 ) ) *traffic = 0;
	return ret_value;
}

/*Function add*/
void  add( int count, double complex *conc, double complex **epsilon, int *traffic )
{
	int i, sw;
	if( *traffic == 1 )
	{
		epsilon[count][1] = Term( count, conc, traffic ) / Term( count - 1, conc, traffic );
		sw = 1;
		for( i = count - 1; i > -1; i-- )
		{
			if( sw > 0 )  epsilon[i][1] = epsilon[i + 1][0] + epsilon[i + 1][1] - epsilon[i][0];
			if( sw < 0 )  epsilon[i][1] = ( epsilon[i + 1][0] * epsilon[i + 1][1] ) / epsilon[i][0];
			sw = -sw;
		}
		for( i = 0; i < count + 1; i++ ) epsilon[i][0] = epsilon[i][1];
	}
}
