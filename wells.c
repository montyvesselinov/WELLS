/* ***************************************************************
 *                       WELLS
 *      Code to compute heads and Drawdown in a well field
 *
 * 						VERSION 1.1
 *
 *       Mishra Phoolendra and Velimir V. Vessilinov (2011)
 * 		   Los Alamos National Laboratory , Los Alamos
 *
 *
 *****************************************************************/
#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <time.h>

#include "design.h"
#include "wells.h"

/* Prototypes of the functions in this file */
void Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd, struct hankel_data *hd, struct laplace_data *ld );
void Save_Problem( struct Problem_Info pi, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd );
void Save_Well_S( char *fileName, struct Well_Data *wd );
void Save_Point_S( char *fileName, struct Point_Data *pd );
void Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
void Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd , struct hankel_data *hd, struct laplace_data *ld );
void CalcSave_Points( char *fileName, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd , struct hankel_data *hd, struct laplace_data *ld );
/*void  Calc_Save_Grid( double time, struct Aquifer_Data ad, struct Well_Data *wd, struct Grid_Data *gd, char *fileName );*/
double Theis( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Papadopulos( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Hantush( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Mishra_conf( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Mishra_unc( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Mishra_leaky_unc( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );

double Theis_OLD( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld );
double Theis_unc( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Hantush_leaky( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Hantush_unc( double x, double y, double z1, double z2,  double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld );
double Ei( double u );
double W( double ra, double rb );
double K0( double x );
double BESSJ0( double X );

int matherr( struct exception *except );
double **double_matrix( int maxCols, int maxRows );
void free_double_matrix( double **matrix, int maxCols );
void free_char_matrix( char **matrix, int maxCols );
char **char_matrix( int maxCols, int maxRows );

/* Prototypes of the functions in other files */
double confined_Step_Lap( struct params *p );
double confined_Step_cumulative_Lap( struct params *p, struct Well_Data *wd, int k , double tf );
double confined_user_def_Lap( struct params *p, struct Well_Data *wd, int k , double tf );
double confined_Linear_Lap( struct params *p, struct Well_Data *wd, int k, double tf );
double confined_Linear_superpose_Lap( struct params *p, struct Well_Data *wd, int k, double tf );
double confined_Linear_convolute_Lap( struct params *p );
double unc_Step_Lap_sat( double y, void *params );
double unc_Linear_Lap_sat( double y, void *params );
double leaky_unc_Step_Lap_sat( double y, void *params );
double leaky_unc_Linear_Lap_sat( double y, void *params );


int main( int argn, char *argv[] )
{
	struct Problem_Info pit;
	struct Aquifer_Data adt;
	struct Well_Data wdt;
	struct Point_Data pdt;
	struct hankel_data hdt;
	struct laplace_data ldt;
	char root[80], *dot;
	clock_t start, end;
	double runTime;
	if( argn < 2 )
	{
		printf( "Usage: wells problem_root\n" );
		printf( "Usage: wells in:problem_file out:well_drawdown out:point_drawdown\n" );
		exit( 1 );
	}
// Start the clock
	start = clock();
	printf( "Argument: %s\n", argv[1] );
	strcpy( root, argv[1] );
	dot = strrchr( root, '.' );
	if( dot != NULL && dot[1] != '/' )
	{
		strcpy( pit.file, argv[1] );
		dot[0] = 0;
	}
	else
		sprintf( pit.file, "%s.wells", argv[1] );
	printf( "Root: %s\n", root );
	printf( "Input file name: %s\n", pit.file );
	Load_Problem( &pit, &adt, &wdt, &pdt, &hdt, &ldt );
	sprintf( pit.file, "%s.wells_debug", root );
	Save_Problem( pit, adt, &wdt, &pdt ); // Save problem to debug
	if( argn < 3 ) sprintf( &pit.file[0], "%s.s_point", root );
	else strcpy( pit.file, argv[2] );
	printf( "\nOutput drawdown file for the points: %s\n", pit.file );
	CalcSave_Points( pit.file, adt, &wdt, &pdt, &hdt, &ldt );
	//Save_Point_S( pit.file, &pdt );
	if( argn == 3 )
	{
		strcpy( pit.file, argv[3] );
		printf( "Output drawdown file for the wells: %s\n", pit.file );
		Calc_Wells( adt, &wdt, &hdt, &ldt );
		Save_Well_S( pit.file, &wdt );
	}
// End clock and compute total Runtime for program
	end = clock();
	runTime = ( end - start ) / ( double ) CLOCKS_PER_SEC ;
	printf( "Computations Finished in %g seconds\n", runTime );
// free up all memories acquired:w
//
	free_char_matrix( wdt.id, wdt.nW );
	free( wdt.x ); free( wdt.y ); free( wdt.r ); free( wdt.h );
	free( wdt.m ); free( wdt.k ); free( wdt.S ); free( wdt.B ); free( wdt.sw );
	free_double_matrix( wdt.Q, wdt.nW ); free_double_matrix( wdt.t, wdt.nW ); free_double_matrix( wdt.s, wdt.nW );
	free_char_matrix( pdt.id, pdt.nP );
	free( pdt.x ); free( pdt.y ); free( pdt.h );
	free_double_matrix( pdt.t, pdt.nP ); free_double_matrix( pdt.s, pdt.nP );
	free_double_matrix( adt.kobs, pdt.nP ); free_double_matrix( adt.Sobs, pdt.nP );
	free_double_matrix( adt.kprod, wdt.nW ); free_double_matrix( adt.Sprod, wdt.nW );
	exit( 0 );
}

void Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd, struct hankel_data *hd, struct laplace_data *ld )
{
	FILE *infile;
	char buf[50], *dot, temp_char[50];
	int  i, j;
	( *ad ).co_aquifer = 0;
	( *ad ).sol_type   = 0;
	if( ( infile = fopen( ( *pi ).file, "r" ) ) == NULL )
	{
		printf( "Problem file could not be read!\n" );
		exit( 1 );
		return;
	}
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ":%[^\n]60s\n", ( *pi ).name );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %s\n", temp_char );
	if( strcasestr( temp_char, "conf" ) )( *ad ).co_aquifer = 1;
	if( strcasestr( temp_char, "uncon" ) )( *ad ).co_aquifer = 2;
	if( strcasestr( temp_char, "leaky-conf" ) )( *ad ).co_aquifer = 3;
	if( strcasestr( temp_char, "leaky-uncon" ) )( *ad ).co_aquifer = 4;
	if( ( *ad ).co_aquifer == 0 ) {printf( "ERROR-Specify correct aquifer Type: confined, unconfined, leaky-confined, leaky-unconfined \n" ); exit( 0 );}
	/*Read Solution Type*/
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %s\n", temp_char );
	if( strcasestr( temp_char, "theis" ) )( *ad ).sol_type = 1;
	if( strcasestr( temp_char, "hantush" ) )( *ad ).sol_type = 2;
	if( strcasestr( temp_char, "papad" ) )( *ad ).sol_type = 3;
	if( strcasestr( temp_char, "mishra" ) )( *ad ).sol_type = 4;
	if( strcasestr( temp_char, "theis-trad" ) )( *ad ).sol_type = 5;
	if( ( *ad ).sol_type == 0 ) {printf( "ERROR-Specify correct Solution Type: Theis,Hantush, Papadopulos,Mishra-Neuman \n" ); exit( 0 );}
	printf( "Solution Type is %i \n", ( *ad ).sol_type );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i %i %lf %lf %lf\n", &( *ad ).co_x_axis, &( *ad ).co_y_axis, &( *ad ).origin_x, &( *ad ).origin_y, &( *ad ).origin_angle );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &( *wd ).nW );
	wd->id = char_matrix( ( *wd ).nW, MAXNAME );
	wd->x = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->y = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->r = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->nQ = ( int * ) malloc( ( *wd ).nW * sizeof( int ) );
	wd->h = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->m = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->k = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->S = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->sy = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->akD = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->acD = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->psiD = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->kD = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->CwD = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->B = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->sw = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->Q = ( double * * ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->t = ( double * * ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->s = ( double * * ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->d = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->l = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->l_unsat = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->pType = ( int * ) malloc( ( *wd ).nW * sizeof( double ) );
	printf( "\nNumber of pumping wells: %d\n", ( *wd ).nW );
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		//printf( "Punmping well #%d: ", i + 1 );
		fscanf( infile, "%[^:]s", buf );
		fscanf( infile, ":%s\n", ( *wd ).id[i] );
		printf( "Pumping well #%i: %s", i + 1, ( *wd ).id[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &( *wd ).x[i], &( *wd ).y[i], &( *wd ).r[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf\n",  &( *wd ).d[i], &( *wd ).l[i] );
		printf( " x=%g y=%g rw=%g d=%g l=%g", ( *wd ).x[i], ( *wd ).y[i], ( *wd ).r[i], ( *wd ).d[i], ( *wd ).l[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).h[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).m[i] ); ( *wd ).m[i] = pow( 10, ( *wd ).m[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).k[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).S[i] );
		/* Read Anisotropic Ratio */
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).kD[i] ); ( *wd ).kD[i] = pow( 10, ( *wd ).kD[i] );
		/* Read wellbore Storage Parameter */
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).CwD[i] ); ( *wd ).CwD[i] = pow( 10, ( *wd ).CwD[i] );
		/* Read The Leakage Parameter */
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).B[i] ); ( *wd ).B[i] = pow( 10, ( *wd ).B[i] );
		/* Read Unsaturated Parameters */
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).l_unsat[i] ); ( *wd ).l_unsat[i] = pow( 10, ( *wd ).l_unsat[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).sy[i] ); ( *wd ).sy[i] = pow( 10, ( *wd ).sy[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).akD[i] ); ( *wd ).akD[i] = pow( 10, ( *wd ).akD[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).acD[i] ); ( *wd ).acD[i] = pow( 10, ( *wd ).acD[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *wd ).psiD[i] );
		/* Read Linear or step Variation */
		//fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %d\n", &(*wd).pType[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %s\n", temp_char );
		if( strcasestr( temp_char, "step" ) )( *wd ).pType[i] = 0;
		if( strcasestr( temp_char, "linear" ) )( *wd ).pType[i] = 1;
		if( strcasestr( temp_char, "linear-super" ) )( *wd ).pType[i] = 2;
		if( strcasestr( temp_char, "linear-conv" ) )( *wd ).pType[i] = 3;
		if( strcasestr( temp_char, "step-lap" ) )( *wd ).pType[i] = 4;
		if( strcasestr( temp_char, "user-def" ) )( *wd ).pType[i] = 5;
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i \n",  &( *wd ).nQ[i] );
		printf( " Q=%i ", ( *wd ).nQ[i] );
		wd->Q[i] = ( double * ) malloc( ( *wd ).nQ[i] * sizeof( double ) );
		wd->t[i] = ( double * ) malloc( ( *wd ).nQ[i] * sizeof( double ) );
		wd->s[i] = ( double * ) malloc( ( *wd ).nQ[i] * sizeof( double ) );
		printf( "h0=%g m=%g k=%g S=%g kD=%g CwD=%g B=%g\n", ( *wd ).h[i], ( *wd ).m[i], pow( 10, ( *wd ).k[i] ), pow( 10, ( *wd ).S[i] ), ( *wd ).kD[i], ( *wd ).CwD[i], ( *wd ).B[i] );
		if( ( *wd ).pType[i] == 0 ) printf( "Using Step varying pumping rates\n" );
		if( ( *wd ).pType[i] == 1 ) printf( "Using Linearly varying pumping rates\n" );
		if( ( *wd ).pType[i] == 4 ) printf( "Using Step changes with Lapalce Inverison\n" );
		if( ( *wd ).pType[i] == 5 ) printf( "Using User Defined Function for Pumping Rate\n" );
		for( j = 0; j < ( *wd ).nQ[i]; j++ )
			fscanf( infile, "%lf %lf\n", &( *wd ).t[i][j], &( *wd ).Q[i][j] );
	}
	/* Read Numerical Inversion parameters */
	fgets( buf, 50, infile );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i %i %lf %lf %i \n", &( *hd ).nZeros_min, &( *hd ).nZeros_max, &( *hd ).EPSABS, &( *hd ).EPSREL, &( *hd ).NUMITER );
	printf( "nZeros_min = %d nZeros_max=%d EPSABS=%g EPSREL = %g NUMITER =%d \n", ( *hd ).nZeros_min, ( *hd ).nZeros_max, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf %lf \n", &( *ld ).DHTOL, &( *ld ).DHALPHA, &( *ld ).time_max );
	/* Point data */
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &( *pd ).nP );
	printf( "\nNumber of observation points: %d\n", ( *pd ).nP );
	pd->id = char_matrix( ( *pd ).nP, MAXNAME );
	pd->x = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	pd->y = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	pd->z1 = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	pd->z2 = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	pd->nT = ( int * ) malloc( ( *pd ).nP * sizeof( int ) );
	pd->t = ( double ** ) malloc( ( *pd ).nP * sizeof( double * ) );
	pd->s = ( double ** ) malloc( ( *pd ).nP * sizeof( double * ) );
	pd->h = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	pd->slope = ( double * ) malloc( ( *pd ).nP * sizeof( double ) );
	// Allocate Number of Zeros for Integration
	pd->nZeros = ( int * ) malloc( ( *pd ).nP * sizeof( double ) );
	for( i = 0; i < ( *pd ).nP; i++ )
	{
		printf( "Observation point #%d: ", i + 1 );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %s \n", ( *pd ).id[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf \n", &( *pd ).x[i], &( *pd ).y[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf \n", &( *pd ).z1[i], &( *pd ).z2[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *pd ).h[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %s\n", buf );
		dot = strstr( buf, "D" );
		if( dot != NULL ) { dot[0] = 'e'; }
		sscanf( buf, "%lf", &( *pd ).slope[i] );
		( *pd ).slope[i] = pow( 10, ( *pd ).slope[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &( *pd ).nT[i] );
		pd->t[i] = ( double * ) malloc( ( *pd ).nT[i] * sizeof( double ) );
		pd->s[i] = ( double * ) malloc( ( *pd ).nT[i] * sizeof( double ) );
		printf( "%s x=%g y=%g z1=%g z2=%g nT=%i h0=%g i=%g\n", ( *pd ).id[i], ( *pd ).x[i], ( *pd ).y[i], ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).nT[i], ( *pd ).h[i], ( *pd ).slope[i] );
		for( j = 0; j < ( *pd ).nT[i]; j++ )
			fscanf( infile, "%lf\n", &( *pd ).t[i][j] );
	}
	fclose( infile );
	ad->kobs = ( double ** ) double_matrix( ( *pd ).nP, ( *wd ).nW );
	ad->Sobs = ( double ** ) double_matrix( ( *pd ).nP, ( *wd ).nW );
	for( i = 0; i < ( *pd ).nP; i++ )
		for( j = 0; j < ( *wd ).nW; j++ )
		{
			( *ad ).kobs[i][j] = pow( 10, ( *wd ).k[j] );
			( *ad ).Sobs[i][j] = pow( 10, ( *wd ).S[j] );
		}
	( *ad ).kprod = double_matrix( ( *wd ).nW, ( *wd ).nW );
	( *ad ).Sprod = double_matrix( ( *wd ).nW, ( *wd ).nW );
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		for( j = 0; j < ( *wd ).nW; j++ )
		{
			( *ad ).kprod[j][i] = pow( 10, ( ( *wd ).k[j] + ( *wd ).k[i] ) / 2 );
			( *ad ).Sprod[j][i] = pow( 10, ( ( *wd ).S[j] + ( *wd ).S[i] ) / 2 );
		}
	}
}

void Save_Problem( struct Problem_Info pi, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd )
{
	FILE *outfile;
	int  i, j;
	if( ( outfile = fopen( pi.file, "w" ) ) == NULL )
	{
		printf( "Problem file could not be opened!\n" );
		return;
	}
	fprintf( outfile, "Problem name: %s\n", pi.name );
	fprintf( outfile, "Aquifer type        : %i\n", ad.co_aquifer );
	/*Print Solution Type*/
	if( ad.sol_type == 1 ) fprintf( outfile, "Solution type        : Theis \n" );
	if( ad.sol_type == 2 ) fprintf( outfile, "Solution type        : Hantush \n" );
	if( ad.sol_type == 3 ) fprintf( outfile, "Solution type        : Papadopulos & Cooper \n" );
	if( ad.sol_type == 4 ) fprintf( outfile, "Solution type        : Mishra-Neuman \n" );
	fprintf( outfile, "Boundary information: %i %i %g %g %g\n", ad.co_x_axis, ad.co_y_axis, ad.origin_x, ad.origin_y, ad.origin_angle );
	fprintf( outfile, "--- Number of  wells: %i\n", ( *wd ).nW );
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		fprintf( outfile, "%s %g %g %g %i\n", ( *wd ).id[i], ( *wd ).x[i], ( *wd ).y[i], ( *wd ).r[i], ( *wd ).nQ[i] );
		fprintf( outfile, "Initial Head        : %g\n", ( *wd ).h[i] );
		fprintf( outfile, "Aquifer Tickness    : %g\n", ( *wd ).m[i] );
		fprintf( outfile, "Permeability        : %g\n", ( *wd ).k[i] );
		fprintf( outfile, "Storage coefficient : %g\n", ( *wd ).S[i] );
		fprintf( outfile, "Leakage coefficient : %g\n", ( *wd ).B[i] );
		for( j = 0; j < ( *wd ).nQ[i]; j++ )
			fprintf( outfile, "%g %g\n", ( *wd ).t[i][j], ( *wd ).Q[i][j] );
	}
	fprintf( outfile, "--- Number of points: %i\n", ( *pd ).nP );
	for( i = 0; i < ( *pd ).nP; i++ )
	{
		fprintf( outfile, "%s %g %g %g %i\n", ( *pd ).id[i], ( *pd ).x[i], ( *pd ).y[i], ( *pd ).h[i], ( *pd ).nT[i] );
		fprintf( outfile, "Initial Head        : %g\n", ( *pd ).h[i] );
		fprintf( outfile, "Water-level slope   : %g\n", ( *pd ).slope[i] );
		/*		fprintf( outfile, "Permeability        : %g\n", (*pd).k[i] );
				fprintf( outfile, "Storage coefficient : %g\n", (*wd).S[i] );*/
		for( j = 0; j < ( *pd ).nT[i]; j++ )
			fprintf( outfile, "%g\n", ( *pd ).t[i][j] );
	}
	fprintf( outfile, "Prod well, obs well, Harm. mean perm., geom. mean storage coeff.\n" );
	for( i = 0; i < ( *pd ).nP; i++ )
		for( j = 0; j < ( *wd ).nW; j++ )
			fprintf( outfile, "%i %i %g %g\n", i, j, ad.kobs[i][j], ad.Sobs[i][j] );
	fprintf( outfile, "Prod well, prod well, Harm. mean perm., geom. mean storage coeff.\n" );
	for( i = 0; i < ( *wd ).nW; i++ )
		for( j = 0; j < ( *wd ).nW; j++ )
			fprintf( outfile, "%i %i %g %g\n", i, j, ad.kprod[i][j], ad.Sprod[i][j] );
	fclose( outfile );
}

void Save_Well_S( char *fileName, struct Well_Data *wd )
{
	FILE *outfile;
	int i, j;
	if( ( outfile = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Well drawdowns could not be saved!\n" );
		return;
	}
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		fprintf( outfile, "%s %g %g %i\n", ( *wd ).id[i], ( *wd ).x[i], ( *wd ).y[i], ( *wd ).nQ[i] );
		for( j = 0; j < ( *wd ).nQ[i]; j++ )
			fprintf( outfile, "%.7g %g %g %.6f\n", ( *wd ).t[i][j], ( *wd ).Q[i][j], ( *wd ).s[i][j], ( *wd ).h[i] - ( *wd ).s[i][j] );
	}
	fclose( outfile );
}

void Save_Point_S( char *fileName, struct Point_Data *pd )
{
	FILE *outfile;
	int i, j;
	if( ( outfile = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Data could not be saved!\n" );
		return;
	}
	for( i = 0; i < ( *pd ).nP; i++ )
	{
		fprintf( outfile, "%s %g %g %i\n", ( *pd ).id[i], ( *pd ).x[i], ( *pd ).y[i], ( *pd ).nT[i] );
		for( j = 0; j < ( *pd ).nT[i]; j++ )
			fprintf( outfile, "%.7g %g %.6f\n", ( *pd ).t[i][j], ( *pd ).s[i][j], ( *pd ).h[i] - ( *pd ).s[i][j] );
	}
	fclose( outfile );
}

void Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	int   i, j, code1, code2, code3;
	double x, y, xc, yc, ca, sa, z, ( *fp )();
	// TODO this function is currently not used
	switch( ad.co_aquifer )
	{
		case CONFINED:
			switch( ad.sol_type )
			{
				case 1:
					fp = Theis; printf( "Using Theis (1935) Solution For Confined Aquifer \n" ); break;
				case 2:
					fp = Hantush; printf( "Using Hantush (1964) Solution For Confined Aquifer \n" ); break;
				case 3:
					fp = Papadopulos; printf( "Using Papadopulos and Cooper (1967) Solution For Confined Aquifer \n" ); break;
				case 4:
					fp = Mishra_conf; printf( "Using Mishra-Neuman (2011) Solution For Confined Aquifer \n" ); break;
				case 5:
					fp = Theis_OLD; printf( "Using Traditional Theis (no-Laplace Inversion)\n" ); break;
			} break;
			//case UNCONFINED:     fp = Theis_unc; break; TODO add this for fun
		case UNCONFINED:
			switch( ad.sol_type )
			{
				case 1:
					fp = Theis_unc; printf( "Using Theis (1935) solution for Unconfined Aquifers\n" ); break;
				case 2:
					fp = Hantush_unc; printf( "Using Hantush (1964) solution for Unconfined Aquifers\n" ); break;
				case 4:
					fp = Mishra_unc; printf( "Using Mishra-Neuman (2011) solution for Unconfined Aquifers\n" ); break;
			}
		case LEAKY:             fp = Hantush_leaky; break;
		case LEAKY_UNCONFINED: // leaky-unconfined
			switch( ad.sol_type )
			{
				case 2:
					fp = Hantush_unc; printf( "Using Hantush solution for unconfined aquifers\n" ); break;
				case 4:
					fp = Mishra_leaky_unc; printf( "Using Mishra et.al. 2011 solution for leaky-unconfined aquifers\n" ); break;
			}
	}
	if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
	{
		ca = cos( ad.origin_angle );
		sa = sin( ad.origin_angle );
		code1 = code2 = code3 = 0;
		if( ad.co_x_axis != NO_BOUNDARY ) code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_y_axis != NO_BOUNDARY ) code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY ) code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	}
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		z = ( *wd ).d[i] + ( ( *wd ).l[i] - ( *wd ).d[i] ) / 2;
		for( j = 0; j < ( *wd ).nQ[i]; j++ )
		{
			( *wd ).s[i][j] = fp( ( *wd ).x[i], ( *wd ).y[i], z, ( *wd ).t[i][j], ( *wd ).h[i], ad.kprod[i], ad.Sprod[i], wd, hd , ld );
			if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
			{
				xc = ( *wd ).x[i] - ad.origin_x;
				yc = ( *wd ).y[i] - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
				if( code1 )
				{
					xc = x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + x * sa + ad.origin_y;
					( *wd ).s[i][j] += code1 * fp( xc, yc, z, ( *wd ).t[i][j], ( *wd ).h[i], ad.kprod[i], ad.Sprod[i], wd, hd , ld );
				}
				if( code2 )
				{
					xc = -x * ca - y * sa + ad.origin_x;
					yc = y * ca + -x * sa + ad.origin_y;
					( *wd ).s[i][j] += code2 * fp( xc, yc, z, ( *wd ).t[i][j], ( *wd ).h[i], ad.kprod[i], ad.Sprod[i], wd , hd, ld );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					( *wd ).s[i][j] += code3 * fp( xc, yc, z, ( *wd ).t[i][j], ( *wd ).h[i], ad.kprod[i], ad.Sprod[i], wd, hd, ld );
				}
			}
		}
	}
}

void Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd , struct hankel_data *hd, struct laplace_data *ld )
{
	int   i, j, code1, code2, code3;
	double x, y, xc, yc, ca, sa, ( *fp )();
	switch( ad.co_aquifer )
	{
		case CONFINED:
			switch( ad.sol_type )
			{
				case 1:
					fp = Theis; printf( "Using Theis(1935) Solution For Confined Aquifer \n" ); break;
				case 2:
					fp = Hantush; printf( "Using Hantush(1964) Solution For Confined Aquifer \n" ); break;
				case 3:
					fp = Papadopulos; printf( "Using Papadopulos and Cooper (1967) Solution For Confined Aquifer \n" ); break;
				case 4:
					fp = Mishra_conf; printf( "Using Mishra-Neuman (2011) Solution For Confined Aquifer \n" ); break;
				case 5:
					fp = Theis_OLD; printf( "Using Traditional Theis (no-Laplace Inversion)\n" ); break;
			} break;
			//case UNCONFINED:     fp = Theis_unc; break;
		case UNCONFINED:
			switch( ad.sol_type )
			{
				case 1:
					fp = Theis_unc; printf( "Using Theis (1935) solution for Unconfined Aquifers\n" ); break;
				case 2:
					fp = Hantush_unc; printf( "Using Hantush (1964) solution for Unconfined Aquifers\n" ); break;
				case 4:
					fp = Mishra_unc; printf( "Using Mishra-Neuman (2011) solution for Unconfined Aquifers\n" ); break;
			}
		case LEAKY:             fp = Hantush_leaky; break;
		case LEAKY_UNCONFINED:
			switch( ad.sol_type )
			{
				case 2:
					fp = Hantush_unc; printf( "Using Hantush solution for unconfined aquifers\n" ); break;
				case 4:
					fp = Mishra_leaky_unc; printf( "Using Mishra et.al. 2011 solution for leaky-unconfined aquifers\n" ); break;
			}
	}
	if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
	{
		ca = cos( ad.origin_angle );
		sa = sin( ad.origin_angle );
		code1 = code2 = code3 = 0;
		if( ad.co_x_axis != NO_BOUNDARY ) code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_y_axis != NO_BOUNDARY ) code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY ) code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	}
	for( i = 0; i < ( *pd ).nP; i++ )
	{
		for( j = 0; j < ( *pd ).nT[i]; j++ )
		{
			( *pd ).s[i][j] = ( ( *pd ).t[i][j] - ( *pd ).t[i][0] ) * ( *pd ).slope[i] + fp( ( *pd ).x[i], ( *pd ).y[i], ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd , hd, ld );
			if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
			{
				xc = ( *pd ).x[i] - ad.origin_x;
				yc = ( *pd ).y[i] - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
				if( code1 )
				{
					xc = x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + x * sa + ad.origin_y;
					( *pd ).s[i][j] += code1 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd , hd, ld );
				}
				if( code2 )
				{
					xc = -x * ca - y * sa + ad.origin_x;
					yc = y * ca + -x * sa + ad.origin_y;
					( *pd ).s[i][j] += code2 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd , hd, ld );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					( *pd ).s[i][j] += code3 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd , hd, ld );
				}
			}
		}
	}
}

void CalcSave_Points( char *fileName, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd , struct hankel_data *hd, struct laplace_data *ld )
{
	int   i, j, k, code1, code2, code3;
	double x, y, xc, yc, ca, sa,
		   ( *fp )( double x, double y, double z1, double z2 , double time, double h, double * k, double * S, struct Well_Data * wd, struct hankel_data * hd, struct laplace_data * ld );
	FILE *outfile;
	if( ( outfile = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Data could not be saved!\n" );
		return;
	}
	switch( ad.co_aquifer )
	{
		case CONFINED:
			switch( ad.sol_type )
			{
				case 1:
					fp = Theis; printf( "Using Theis(1935) Solution For Confined Aquifer \n" ); break;
				case 2:
					fp = Hantush; printf( "Using Hantush(1964) Solution For Confined Aquifer \n" ); break;
				case 3:
					fp = Papadopulos; printf( "Using Papadopulos and Cooper (1967) Solution For Confined Aquifer \n" ); break;
				case 4:
					fp = Mishra_conf; printf( "Using Mishra-Neuman (2011) Solution For Confined Aquifer \n" ); break;
				case 5:
					fp = Theis_OLD; printf( "Using Traditional Theis (no-Laplace Inversion)\n" ); break;
			} break;
			//case UNCONFINED:       fp = Theis_unc; break; TODO we need to be able to envoke this for simple test analyses 
		case UNCONFINED:        fp = Mishra_unc; printf( "Using Mishra and Neuman (2011) Solution for Unconfined aquifers\n" ); break;
		case LEAKY:             fp = Hantush; printf( "Using Hantush (1964) Solution For Leaky Confined Aquifer \n" ); break;
		case LEAKY_UNCONFINED:
			switch( ad.sol_type )
			{
				case 4:
					fp = Mishra_leaky_unc; printf( "Using Mishra et.al. 2011 solution for leaky-unconfined aquifers\n" ); break;
			}
	}
	if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
	{
		ca = cos( ad.origin_angle );
		sa = sin( ad.origin_angle );
		code1 = code2 = code3 = 0;
		if( ad.co_x_axis != NO_BOUNDARY ) code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_y_axis != NO_BOUNDARY ) code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY ) code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	}
	for( i = 0; i < ( *wd ).nW; i++ )( *wd ).sw[i] = 0;
	for( i = 0; i < ( *pd ).nP; i++ )
	{
		fprintf( outfile, "[%s x:%g y:%g #times:%i] Cols: Time elev total_dd trend_dd\n", ( *pd ).id[i], ( *pd ).x[i], ( *pd ).y[i], ( *pd ).nT[i] );
		printf( "[%s x:%g y:%g #times:%i] Cols: Time elev total_dd trend_dd\n", ( *pd ).id[i], ( *pd ).x[i], ( *pd ).y[i], ( *pd ).nT[i] );
		for( j = 0; j < ( *pd ).nT[i]; j++ )
		{
			( *pd ).s[i][j] = ( ( *pd ).t[i][j] - ( *pd ).t[i][0] ) * ( *pd ).slope[i] + fp( ( *pd ).x[i], ( *pd ).y[i], ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd , hd, ld );
			if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
			{
				xc = ( *pd ).x[i] - ad.origin_x;
				yc = ( *pd ).y[i] - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
				if( code1 )
				{
					xc = x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + x * sa + ad.origin_y;
					( *pd ).s[i][j] += code1 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd, hd, ld );
				}
				if( code2 )
				{
					xc = -x * ca - y * sa + ad.origin_x;
					yc = y * ca + -x * sa + ad.origin_y;
					( *pd ).s[i][j] += code2 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd, hd, ld );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					( *pd ).s[i][j] += code3 * fp( xc, yc, ( *pd ).z1[i], ( *pd ).z2[i], ( *pd ).t[i][j], ( *pd ).h[i], ad.kobs[i], ad.Sobs[i], wd, hd, ld );
				}
			}
			fprintf( outfile, "%g %.6f %g %g", ( *pd ).t[i][j], ( *pd ).h[i] - ( *pd ).s[i][j], ( *pd ).s[i][j], ( ( *pd ).t[i][j] - ( *pd ).t[i][0] ) * ( *pd ).slope[i] );
			printf( "%g %.6f %g %g", ( *pd ).t[i][j], ( *pd ).h[i] - ( *pd ).s[i][j], ( *pd ).s[i][j], ( ( *pd ).t[i][j] - ( *pd ).t[i][0] ) * ( *pd ).slope[i] );
			for( k = 0; k < ( *wd ).nW; k++ )
			{
				fprintf( outfile, " %g", ( *wd ).sw[k] );
				printf( " %g", ( *wd ).sw[k] );
				( *wd ).sw[k] = 0;
			}
			fprintf( outfile, "\n" );
			printf( "\n" );
			fflush( outfile );
		}
	}
	fclose( outfile );
}
/*   BELOW ARE INDIVIDUAL SOLUTIONS */

double Mishra_conf( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, tf, t0, td, T;
	int i, n, nR;
	struct params p;
	s = 0.0; td = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		p.rwD = sqrt( rw / r2 );
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		p.lD = ( *wd ).l[i] / ( *wd ).m[i];
		p.dD = ( *wd ).d[i] / ( *wd ).m[i];
		p.zD1 = z1 / ( *wd ).m[i]; p.zD2 = z2 / ( *wd ).m[i];
		p.cwD = ( *wd ).CwD[i];
		p.kD = ( *wd ).kD[i];
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				s += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s;
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					s += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					( *wd ).sw[i] = s;
				}
			}
			// printf(" %g %g \n",time, s);
		}
		if( ( *wd ).pType[i] == 4 )
		{
			p.ts = time * tf;
			s = confined_Step_cumulative_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			p.ts = time * tf;
			s = confined_Linear_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
			// printf(" %g  %g \n",time, s);
		} // End Linear Case
		if( ( *wd ).pType[i] == 2 )
		{
			p.ts = time * tf;
			s = confined_Linear_superpose_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case with superposition
		if( ( *wd ).pType[i] == 3 )
		{
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) // loop through all the pumping rates
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				s += confined_Linear_convolute_Lap( &p ) / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s ;
			} // End Loop over observation times
			printf( " %g  %g \n", time, s );
		} // End Linear Case
	} // End Step or Linear Condition
	return( s );
}

double Theis( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, tf, t0, td, sub1, T;
	int i, n, nR;
	struct params p;
	s = 0.0; td = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		p.rwD = 0.0;
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		sub1 = ( double ) 4.0 * M_PI * T;
		p.lD = 1.0;
		p.dD = 0.0;
		p.zD1 = 0.5; p.zD2 = 0.5;
		p.cwD = 0.0;
		p.kD = 1.0;
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				s += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / sub1;
				( *wd ).sw[i] = s;
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					s += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / sub1;
					( *wd ).sw[i] = s;
				}
			}
			//printf(" %g %g \n",time, s);
		}
		if( ( *wd ).pType[i] == 5 )
		{
			p.ts = ( time - ( *wd ).t[i][0] ) * tf;
			s = confined_user_def_Lap( &p, wd, i, tf ) / sub1;
			( *wd ).sw[i] = s ;
			printf( " %g %g \n", time, s );
		}
		if( ( *wd ).pType[i] == 4 )
		{
			p.ts = ( time - ( *wd ).t[i][0] ) * tf;
			//p.ts = time*tf;
			s = confined_Step_cumulative_Lap( &p, wd, i, tf ) / sub1;
			( *wd ).sw[i] = s ;
			//printf(" %g %g \n",time, s);
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			p.ts = time * tf;
			//p.ts = (time-(*wd).t[i][0])*tf;
			//s = confined_Step_Lap1(&p, wd, i, tf)/sub1; //((double)4.0*M_PI*T);
			s = confined_Linear_Lap( &p, wd, i, tf ) / sub1; //((double)4.0*M_PI*T);
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case -direct method
		if( ( *wd ).pType[i] == 2 )
		{
			p.ts = time * tf;
			s = confined_Linear_superpose_Lap( &p, wd, i, tf ) / sub1; //((double)4.0*M_PI*T);
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case with superposition
		if( ( *wd ).pType[i] == 3 )
		{
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) // loop through all the pumping rates
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				s += confined_Linear_convolute_Lap( &p ) / sub1; //((double)4.0*M_PI*T);
				( *wd ).sw[i] = s ;
			} // End Loop over observation times
			//printf(" %g  %g \n",time, s);
		} // End Linear Case with convolution
	} // End Step or Linear Condition
	return( s );
}

double Hantush( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, tf, t0, td, T, sub1;
	int i, n, nR;
	struct params p;
	s = 0.0; td = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		sub1 = ( double )4.0 * M_PI * T;
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		p.rwD = 0.0 ;
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		p.lD = ( *wd ).l[i] / ( *wd ).m[i];
		p.dD = ( *wd ).d[i] / ( *wd ).m[i];
		p.zD1 = z1 / ( *wd ).m[i]; p.zD2 = z2 / ( *wd ).m[i];
		p.cwD = ( *wd ).CwD[i];
		p.kD = ( *wd ).kD[i];
		//printf(" %g %g \n",p.cwD, p.kD);
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				s += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / sub1;
				( *wd ).sw[i] = s;
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					s += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / sub1;
					( *wd ).sw[i] = s;
				}
			}
			// printf(" %g %g \n",time, s);
		}
		if( ( *wd ).pType[i] == 4 )
		{
			p.ts = time * tf;
			s = confined_Step_cumulative_Lap( &p, wd, i, tf ) / sub1;
			( *wd ).sw[i] = s ;
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			p.ts = time * tf;
			s = confined_Linear_Lap( &p, wd, i, tf ) / sub1;
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case
		if( ( *wd ).pType[i] == 2 )
		{
			p.ts = time * tf;
			s = confined_Linear_superpose_Lap( &p, wd, i, tf ) / sub1;
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case with superposition
		if( ( *wd ).pType[i] == 3 )
		{
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) // loop through all the pumping rates
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				s += confined_Linear_convolute_Lap( &p ) / sub1;
				( *wd ).sw[i] = s ;
			} // End Loop over observation times
			printf( " %g  %g \n", time, s );
		} // End Linear Case with convolutions
	} // End Step or Linear Condition
	return( s );
}

double Papadopulos( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, tf, t0, td, T;
	int i, n, nR;
	struct params p;
	s = 0.0; td = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		p.rwD = sqrt( rw / r2 );
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		p.lD = 1.0 ;
		p.dD = 0.0 ;
		p.zD1 = z1 / ( *wd ).m[i]; p.zD2 = z1 / ( *wd ).m[i];
		p.cwD = ( *wd ).CwD[i];
		p.kD = 1.0;
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				s += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s;
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					s += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					( *wd ).sw[i] = s;
				}
			}
			// printf(" %g %g \n",time, s);
		}
		if( ( *wd ).pType[i] == 4 )
		{
			p.ts = time * tf;
			s = confined_Step_cumulative_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			p.ts = time * tf;
			s = confined_Linear_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case
		if( ( *wd ).pType[i] == 2 )
		{
			p.ts = time * tf;
			s = confined_Linear_superpose_Lap( &p, wd, i, tf ) / ( ( double )4.0 * M_PI * T );
			( *wd ).sw[i] = s ;
			//printf(" %g  %g \n",time, s);
		} // End Linear Case with superposition
		if( ( *wd ).pType[i] == 3 )
		{
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) // loop through all the pumping rates
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				s += confined_Linear_convolute_Lap( &p ) / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s ;
			} // End Loop over observation times
			printf( " %g  %g \n", time, s );
		} // End Linear Case with convolution
	} // End Step or Linear Condition
	return( s );
}

double Mishra_unc( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, sC, sU, tf, sd, t0, td, T;
	int i, n, nR, status, nZeros, jerr[10001];
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( ( *hd ).NUMITER ); 
	gsl_function F;
	double result, error, j0zero[1001];
	void ROOTJ( int N, int NK, double * JZERO, int * IER );
	struct params p;
	sC = 0.0; sU = 0.0; tf = 0.0; sd = 0.0; td = 0.0; s = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		// printf("T=%g S=%g Sy=%g akD=%g \n",T,S[i],(*wd).sy[i],(*wd).akD[i]);
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		p.lD = ( *wd ).l[i] / ( *wd ).m[i];
		p.dD = ( *wd ).d[i] / ( *wd ).m[i];
		p.zD1 = z1 / ( *wd ).m[i]; p.zD2 = z2 / ( *wd ).m[i];
		p.cwD = ( *wd ).CwD[i]; p.rwD = sqrt( rw / r2 );
		p.akD = ( *wd ).akD[i]; p.acD = ( *wd ).acD[i]; p.kD = ( *wd ).kD[i];
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		p.sigma = S[i] / ( *wd ).sy[i]; p.psiD = ( *wd ).psiD[i];
		p.l_unsat = ( *wd ).l_unsat[i];
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
		//     printf("rwD=%g cwD=%g \n",p.rwD,p.cwD);
		/* Compute the Zeros of Bessel functions J0 */
		//nZeros=(int)(25.0/(sqrt(p.kD)*p.rD));
		nZeros = 5;
		ROOTJ( 0, nZeros, j0zero, jerr );
		j0zero[0] = 0.0;
		//printf(" NUMBER OF ZEROS %g \n", j0zero[5]/(sqrt(p.kD)*p.rD));
		//for (n=0;n< nZeros+1;n++) {printf("%d %g \n",n, j0zero[n]/(sqrt(p.kD)*p.rD));j0zero[n]=j0zero[n]/(sqrt(p.kD)*p.rD);}
		for( n = 0; n < nZeros + 1; n++ ) j0zero[n] = j0zero[n] / ( sqrt( p.beta2 ) );
		/* Hankel inversion on Laplace inversed function */
		gsl_set_error_handler_off();
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			F.params = &p;
			F.function = &unc_Step_Lap_sat;
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				sC += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				do
				{
					status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
					//printf("%g %g \n",(*hd).EPSABS, (*wd).m[i]);
					if( nZeros > 1 && isnan( result ) ) nZeros--;
				}
				while( isnan( result ) );
				sU += result * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s = sC + sU;
				//printf("sU=%g \n",sU);
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					sC += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					do
					{
						status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
						//printf("%g \n",status);
						if( ( nZeros > 1 ) && ( isnan( result ) ) ) nZeros--;
					}
					while( isnan( result ) );
					sU += result * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					( *wd ).sw[i] = s = sC + sU;
				}
			}
			//printf(" %g %g %g %g \n",time,sC,sU, s);
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			F.params = &p;
			F.function = &unc_Linear_Lap_sat;
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) /* loop through all the pumping rates */
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				sC += confined_Linear_convolute_Lap( &p ) / ( ( double )4.0 * M_PI * T );
				do
				{
					status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
					if( nZeros > 1 && isnan( result ) ) nZeros = nZeros - 1;
				}
				while( isnan( result ) );
				sU += result / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s = sC + sU;
			} // End Loop over observation times
			printf( " %g %g %g %g \n", time, sC, sU, s );
		} // End Linear Case
	} // End Step or Linear Condition
	gsl_integration_workspace_free( w );
	return( s );
}
/***************************************************/
/*   Leaky Unconfined Aquifer Mishra et.al. 2011   */
/***************************************************/

double Mishra_leaky_unc( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd , struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r2, rw, s, sC, sU, tf, sd, t0, td, T;
	int i, n, nR, status, nZeros, jerr[10001];
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( ( *hd ).NUMITER ); gsl_function F;
	double result, error, j0zero[1001];
	void ROOTJ( int N, int NK, double * JZERO, int * IER );
	struct params p;
	sC = 0.0; sU = 0.0; tf = 0.0; sd = 0.0; td = 0.0; s = 0.0;
	/*Note sd=4*PI*T*s similar to dimesnionless drawdown*/
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		//printf("T=%g S=%g Sy=%g akD=%g \n",T,S[i],(*wd).sy[i],(*wd).akD[i]);
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r2 = dx * dx + dy * dy; /* squared radial distance */
		rw = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r2 < rw ) r2 = rw;  /* if within pumping well compute at well face*/
		tf = T / ( S[i] * r2 ); /* Dimensionless Time factor based on td=(Kr/Ss)t/r2  where S is storativity*/
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		p.lD = ( *wd ).l[i] / ( *wd ).m[i];
		p.dD = ( *wd ).d[i] / ( *wd ).m[i];
		p.zD1 = z1 / ( *wd ).m[i]; p.zD2 = z2 / ( *wd ).m[i];
		p.cwD = ( *wd ).CwD[i]; p.rwD = sqrt( rw / r2 );
		p.akD = ( *wd ).akD[i]; p.acD = ( *wd ).acD[i]; p.kD = ( *wd ).kD[i];
		p.beta2 = ( *wd ).kD[i] * r2 / ( ( *wd ).m[i] * ( *wd ).m[i] );
		p.sigma = S[i] / ( *wd ).sy[i]; p.psiD = ( *wd ).psiD[i];
		p.l_unsat = ( *wd ).l_unsat[i];
		p.DHTOL = ( *ld ).DHTOL; p.DHALPHA = ( *ld ).DHALPHA; p.time_max = ( *ld ).time_max;
		//     printf("rwD=%g cwD=%g \n",p.rwD,p.cwD);
		/* Compute the Zeros of Bessel functions J0 */
		//nZeros=(int)(25.0/(sqrt(p.kD)*p.rD));
		nZeros = 5;
		ROOTJ( 0, nZeros, j0zero, jerr );
		j0zero[0] = 0.0;
		//printf(" NUMBER OF ZEROS %g \n", j0zero[5]/(sqrt(p.kD)*p.rD));
		//for (n=0;n< nZeros+1;n++) {printf("%d %g \n",n, j0zero[n]/(sqrt(p.kD)*p.rD));j0zero[n]=j0zero[n]/(sqrt(p.kD)*p.rD);}
		for( n = 0; n < nZeros + 1; n++ ) j0zero[n] = j0zero[n] / ( sqrt( p.beta2 ) );
		/* Hankel inversion on Laplace inversed function */
		gsl_set_error_handler_off();
// Step Changes
		if( ( *wd ).pType[i] == 0 )
		{
			t0 = time - ( *wd ).t[i][0];
			F.params = &p;
			F.function = &leaky_unc_Step_Lap_sat;
			if( t0 > 0 )
			{
				p.ts = t0 * tf;
				sC += confined_Step_Lap( &p ) * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				do
				{
					status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
					//printf("%g %g \n",(*hd).EPSABS, (*wd).m[i]);
					if( nZeros > 1 && isnan( result ) ) nZeros--;
				}
				while( isnan( result ) );
				sU += result * ( *wd ).Q[i][0] / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s = sC + sU;
				//printf("sU=%g \n",sU);
			}
			for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
			{
				td = time - ( *wd ).t[i][n];
				if( td > 0 )
				{
					p.ts = td * tf;
					sC += confined_Step_Lap( &p ) * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					do
					{
						status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
						//printf("%g \n",status);
						if( ( nZeros > 1 ) && ( isnan( result ) ) ) nZeros--;
					}
					while( isnan( result ) );
					sU += result * ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) / ( ( double )4.0 * M_PI * T );
					( *wd ).sw[i] = s = sC + sU;
				}
			}
			//printf(" %g %g %g %g \n",time,sC,sU, s);
		}
// Linear Changes
		if( ( *wd ).pType[i] == 1 )
		{
			F.params = &p;
			F.function = &leaky_unc_Linear_Lap_sat;
			for( n = 1; n < nR && ( *wd ).t[i][n - 1] < time; n++ ) /* loop through all the pumping rates */
			{
				p.ts = time * tf;
				p.Q1 = ( *wd ).Q[i][n - 1]; p.Q2 = ( *wd ).Q[i][n];
				p.ts1 = ( *wd ).t[i][n - 1] * tf; p.ts2 = ( *wd ).t[i][n] * tf;
				sC += confined_Linear_convolute_Lap( &p ) / ( ( double )4.0 * M_PI * T );
				do
				{
					status = gsl_integration_qagp( &F, j0zero, nZeros, ( *hd ).EPSABS, ( *hd ).EPSREL, ( *hd ).NUMITER, w, &result, &error );
					if( nZeros > 1 && isnan( result ) ) nZeros--;
				}
				while( isnan( result ) );
				sU += result / ( ( double )4.0 * M_PI * T );
				( *wd ).sw[i] = s = sC + sU;
			} // End Loop over observation times
			printf( " %g %g %g %g \n", time, sC, sU, s );
		} // End Linear Case
	} // End Step or Linear Condition
	gsl_integration_workspace_free( w );
	return( s );
}

/*********************************/
/*********************************/
double Theis_OLD( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r, T, sub1, sub2, sub3, sub4, s, ss;
	int i, n, nR;
	s = 0;
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		( *wd ).sw[i] = 0;
		T = ( *wd ).m[i] * k[i]; /* transmissivity */
		/* printf( "%i %g %g\n", i, k[i], S[i] );  */
		sub1 = 12.5663706 * T; /* 4 * PI = 12.5663706 */
		sub2 = S[i] / ( 4.0 * T );
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r = dx * dx + dy * dy; /* squared radial distance */
		sub3 = ( *wd ).r[i] * ( *wd ).r[i]; /* squared well radius */
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		sub4 = time - ( *wd ).t[i][0];
		if( sub4 > 0 )
		{
			ss = ( *wd ).Q[i][0] * Ei( sub3 / sub4 ) / sub1; /* initial drawdown */
			( *wd ).sw[i] += ss;
			s += ss;
		}
		nR = ( *wd ).nQ[i]; /* number of step pumping rate changes */
		for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ ) /* loop through all the pumping rates */
		{
			sub4 = time - ( *wd ).t[i][n];
			if( sub4 > 0 )
			{
				ss = ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) * Ei( sub3 / sub4 ) / sub1;
				( *wd ).sw[i] += ss;
				s += ss;
			}
			//printf(" drawdown = %g \n",s);
		}
	}
	return( s );
}

double Theis_unc( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r, h2, sub1, sub2, sub3, sub4, s;
	int i, n, nR;
	h2 = h * h;
	s = 0;
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		sub1 = 6.2831853 * k[i];
		sub2 = S[i] / ( 4.0 * k[i] * ( *wd ).m[i] );
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r = dx * dx + dy * dy;
		sub3 = ( *wd ).r[i] * ( *wd ).r[i];
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		sub4 = time - ( *wd ).t[i][0];
		if( sub4 > 0 ) s += h - sqrt( h2 - ( *wd ).Q[i][0] * Ei( sub3 / sub4 ) / sub1 );
		nR = ( *wd ).nQ[i];
		for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
		{
			sub4 = time - ( *wd ).t[i][n];
			if( sub4 > 0 ) s += h - sqrt( h2 - ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) * Ei( sub3 / sub4 ) / sub1 );
		}
	}
	return( s );
}

double Hantush_leaky( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r, T, sub1, sub2, sub3, sub4, sub5, s;
	int i, n, nR;
	s = 0;
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		T = ( *wd ).m[i] * k[i];
		sub1 = 12.5663706 * T;
		sub2 = S[i] / ( 4.0 * T );
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < ( *wd ).r[i] ) r = ( *wd ).r[i];
		sub3 = r * r * sub2;
		sub4 = r / ( *wd ).B[i];
		sub5 = time - ( *wd ).t[i][0];
		if( sub5 > 0 ) s += ( *wd ).Q[i][0] * W( sub3 / sub5, sub4 ) / sub1;
		nR = ( *wd ).nQ[i];
		for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
		{
			sub5 = time - ( *wd ).t[i][n];
			if( sub5 > 0 ) s += ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) * W( sub3 / sub5, sub4 ) / sub1;
		}
	}
	return( s );
}

double Hantush_unc( double x, double y, double z1, double z2, double time, double h, double *k, double *S, struct Well_Data *wd, struct hankel_data *hd, struct laplace_data *ld )
{
	double dx, dy, r, h2, sub1, sub2, sub3, sub4, sub5, s;
	int i, n, nR;
	h2 = h * h;
	s = 0;
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		sub1 = 6.2831853 * k[i];
		sub2 = S[i] / ( 4.0 * k[i] * ( *wd ).m[i] );
		dx = x - ( *wd ).x[i];	dy = y - ( *wd ).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < ( *wd ).r[i] ) r = ( *wd ).r[i];
		sub3 = r * r * sub2;
		sub4 = r / ( *wd ).B[i];
		sub5 = time - ( *wd ).t[i][0];
		if( sub5 > 0 ) s += h - sqrt( h2 - ( *wd ).Q[i][0] * W( ( sub3 / sub5 ), sub4 ) / sub1 );
		nR = ( *wd ).nQ[i];
		for( n = 1; n < nR && ( *wd ).t[i][n] < time; n++ )
		{
			sub5 = time - ( *wd ).t[i][n];
			if( sub5 > 0 ) s += h - sqrt( h2 - ( ( *wd ).Q[i][n] - ( *wd ).Q[i][n - 1] ) * W( ( sub3 / sub5 ), sub4 ) / sub1 );
		}
	}
	return( s );
}

/*  The Theis Function - Ei(x)
 */
double Ei( double u )
{
	double s1, s2;
	if( u > 60.0 || u <= 0.0 ) return( 0.0 );
	if( u < 1.0 )
	{
		s1 = ( ( ( ( 0.00107857 * u - 0.00976004 ) * u + 0.05519968 ) * u - 0.24991055 ) * u + 0.99999193 ) * u - 0.57721566;
		return( s1 - log( u ) );
	}
	else
	{
		s1 = ( ( ( u + 8.57332875 ) * u + 18.05901697 ) * u + 8.63476088 ) * u + 0.267773734;
		s2 = ( ( ( u + 9.57332875 ) * u + 25.63295615 ) * u + 21.09965308 ) * u + 3.958496923;
		return( ( 1 / u ) * ( s1 / s2 ) * exp( -u ) );
	}
}

/*  The Hantush Function
 *  The function is calculated by logarithmic integration using Simpson's rule
 *  source: W.Kinzelbach - Groundwater modelling
 */
double W( double ra, double rb )
{
	double ug, og, hi, sub1, sub2, x2, x4, s2, s4;
	int   i;
	if( ra < 1e-20 )
		return( 2 * K0( rb ) );
	if( ra > 2000.0 )
		return( 0.0 );
	ug = log( ra );
	og = 10.0;
	hi = ( og - ug ) / 24;
	sub1 = rb * rb / 4.0;   sub2 = hi * 2;
	x4 = ug + hi;           x2 = ug;
	s4 = s2 = 0.0;          i = 0;
	while( 1 )
	{
		s4 += exp( - exp( x4 ) - sub1 / exp( x4 ) );
		x4 += sub2;
		if( i == 10 ) break;
		x2 += sub2;
		s2 += exp( - exp( x2 ) - sub1 / exp( x2 ) );
		i++;
	}
	return( hi * ( exp( - exp( ug ) - sub1 / exp( ug ) ) + 4 * s4 + 2 * s2 + exp( - exp( og ) - sub1 / exp( og ) ) ) / 3 );
}

/*  The Bessel Function Ko(x)
 */
double K0( double x )
{
	double x1, x2;
	double s1, s2;
	if( x >= 60.0 )
		return( 0.0 );
	else if( x >= 2.0 )
	{
		x1 = 2.0 / x;
		s1 = ( ( ( ( ( 0.00053208 * x1 - 0.00251540 ) * x1 + 0.00587872 ) * x1 - 0.01062446 ) * x1 + 0.02189568 ) * x1 - 0.07832358 ) * x1 + 1.25331414;
		return( s1 * exp( -x ) / sqrt( x ) );
	}
	else
	{
		x1 = x * x / 4.0;
		x2 = x * x / 14.0625;
		s1 = ( ( ( ( ( 0.00000740 * x1 + 0.00010750 ) * x1 + 0.00262698 ) * x1 + 0.03488590 ) * x1 + 0.23069756 ) * x1 + 0.42278420 ) * x1 - 0.57721566;
		s2 = ( ( ( ( ( 0.0045813 * x2 + 0.0360768 ) * x2 + 0.2659732 ) * x2 + 1.2067492 ) * x2 + 3.0899424 ) * x2 + 3.5156229 ) * x2 + 1.0;
		return( -log( 0.5 * x ) * s2 + s1 );
	}
}

int matherr( struct exception *except )
{
	if( except->type == DOMAIN )
	{
		if( strcmp( except->name, "sqrt" ) == 0 )
		{
			except->retval = 0.0;
			printf( "MATH ERROR: Negative root!" );
			return 1;
		}
		return 0;
	}
	else
	{
		printf( "ERROR!" );
		return 0;
	}
}

char **char_matrix( int maxCols, int maxRows )
{
	char **matrix;
	int i;
	if( ( matrix = ( char ** ) malloc( maxCols * sizeof( char * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if( ( matrix[i] = ( char * ) malloc( maxRows * sizeof( char ) ) ) == NULL )
		{
			for( i--; i >= 0; i-- )
				free( matrix[i] );
			free( matrix );
			return( NULL );
		}
	return( matrix );
}

double **double_matrix( int maxCols, int maxRows )
{
	double **matrix;
	int i;
	if( ( matrix = ( double ** ) malloc( maxCols * sizeof( double * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if( ( matrix[i] = ( double * ) malloc( maxRows * sizeof( double ) ) ) == NULL )
		{
			for( i--; i >= 0; i-- )
				free( matrix[i] );
			free( matrix );
			return( NULL );
		}
	return( matrix );
}

void free_double_matrix( double **matrix, int maxCols )
{
	int i;
	for( i = 0; i < maxCols; i++ )
		free( matrix[i] );
	free( matrix );
}

void free_char_matrix( char **matrix, int maxCols )
{
	int i;
	for( i = 0; i < maxCols; i++ )
		free( matrix[i] );
	free( matrix );
}

/*****************FUNTIONS TO FIND ROOT OF J-Beesel Functions *************/
/* Function to change sign of number required by Root finidng functions */
double Sign( double a, double b )
{
	if( b >= 0.0 )
		return ( fabs( a ) );
	else
		return ( -fabs( a ) );
}

//----------------------------------------------------------------------
void ROOTJ( int N, int NK, double *JZERO, int *IER )
{
	/*----------------------------------------------------------------------
	!     CALCULATE THE FIRST NK ZEROES OF BESSEL FUNCTION J(N,X)
	!
	!     INPUTS:
	!       N    ORDER OF FUNCTION J (INTEGER >= 0)                  I*4
	!       NK   NUMBER OF FIRST ZEROES  (INTEGER > 0)               I*4
	!     OUTPUTS:
	!       JZERO(NK)  TABLE OF FIRST ZEROES (ABCISSAS)              R*8
	!       IER(NK)    TABLE OF ERROR CODES (MUST BE ZEROES)         I*4
	!
	!     REFERENCE :
	!     ABRAMOWITZ M. & STEGUN IRENE A.
	!     HANDBOOK OF MATHEMATICAL FUNCTIONS
	! -------------------------------------------------------------------- */
	double ZEROJ, B0, B1, B2, B3, B5, B7, T0, T1, T3, T5, T7, FN, FK;
	double C1, C2, C3, C4, C5, F1, F2, F3, TOL, ERRJ, PI;
	int IERROR, K, NITMX;
	void SECANT( int N, int NITMX, double TOL, double * ZEROJ, int * IER );
	double BESSJ( int N, double X );
	double BESSJ0( double X );
	double BESSJ1( double X );
	TOL = 1E-8; NITMX = 10; PI = 4.0 * atan( 1.0 );
	C1 = 1.8557571; C2 = 1.033150; C3 = 0.00397; C4 = 0.0908; C5 = 0.043;
	FN = 1.0 * N;
//    FIRST ZERO
	if( N == 0 )
	{
		ZEROJ = C1 + C2 - C3 - C4 + C5;
		SECANT( N, NITMX, TOL, &ZEROJ, &IERROR );
		IER[1] = IERROR;
		JZERO[1] = ZEROJ;
	}
	else
	{
		F1 = pow( FN, ( 1.0 / 3.0 ) );
		F2 = F1 * F1 * FN;
		F3 = F1 * FN * FN;
		ZEROJ = FN + C1 * F1 + ( C2 / F1 ) - ( C3 / FN ) - ( C4 / F2 ) + ( C5 / F3 );
		SECANT( N, NITMX, TOL, &ZEROJ, &IERROR );
		IER[1] = IERROR;
		JZERO[1] = ZEROJ;
	}
	T0 = 4.0 * FN * FN;
	T1 = T0 - 1.0;
	T3 = 4.0 * T1 * ( 7.0 * T0 - 31.0 );
	T5 = 32.0 * T1 * ( ( 83.0 * T0 - 982.0 ) * T0 + 3779.0 );
	T7 = 64.0 * T1 * ( ( ( 6949.0 * T0 - 153855.0 ) * T0 + 1585743.0 ) * T0 - 6277237.0 );
//    OTHER ZEROES
	for( K = 2; K <= NK; K++ )
	{
		JZERO[K] = 0.0;
		FK = 1.0 * K;
//      MAC MAHON'S SERIES FOR K>>N
		B0 = ( FK + 0.5 * FN - 0.25 ) * PI;
		B1 = 8.0 * B0;
		B2 = B1 * B1;
		B3 = 3.0 * B1 * B2;
		B5 = 5.0 * B3 * B2;
		B7 = 7.0 * B5 * B2;
		ZEROJ = B0 - ( T1 / B1 ) - ( T3 / B3 ) - ( T5 / B5 ) - ( T7 / B7 );
		ERRJ = fabs( BESSJ( N, ZEROJ ) );
//      IMPROVE SOLUTION USING PROCEDURE SECANT
		if( ERRJ > TOL ) SECANT( N, NITMX, TOL, &ZEROJ, &IERROR );
		JZERO[K] = ZEROJ;
		IER[K] = IERROR;
	}
}
// ------------------------------------------------------------------------------
void SECANT( int N, int NITMX, double TOL, double *ZEROJ, int *IER )
{
	//Labels: e5,e10,e15,e20
	double P0, P1, Q0, Q1, DP, P;
	double C[3];
	int IT, NEV, NTRY;
	double BESSJ( int N, double X );
	double BESSJ0( double X );
	double BESSJ1( double X );
	C[1] = 0.95; C[2] = 0.999;
	NTRY = 1; P = 0.0;
	*IER = 0;
e5:   P0 = C[NTRY] * ( *ZEROJ );
	P1 = *ZEROJ;
	NEV = 2;
	Q0 = BESSJ( N, P0 );
	Q1 = BESSJ( N, P1 );
	for( IT = 1; IT <= NITMX; IT++ )
	{
		if( Q1 == Q0 ) goto e15;
		P = P1 - Q1 * ( P1 - P0 ) / ( Q1 - Q0 );
		DP = P - P1;
		if( IT == 1 ) goto e10;
		if( fabs( DP ) < TOL ) goto e20;
e10:    NEV = NEV + 1;
		P0 = P1;
		Q0 = Q1;
		P1 = P;
		Q1 = BESSJ( N, P1 );
	}
e15:  NTRY++;
	if( NTRY <= 2 ) goto e5;
	*IER = NTRY;
e20:  *ZEROJ = P;
}

// --------------------------------------------------------------------
double BESSJ( int N, double X )
{
	/*    THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
	      OF ORDER N, INTEGER FOR ANY REAL X. WE USE HERE THE CLASSICAL
	      RECURRENT FORMULA, WHEN  X > N. FOR X < N, THE MILLER'S ALGORITHM
	      IS USED TO AVOID OVERFLOWS.
	      REFERENCE :
	      C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
	      MATHEMATICAL TABLES, VOL.5, 1962.  */
	int IACC = 40;
	double BIGNO = 1e10, BIGNI = 1E-10;
	double TOX, BJM, BJ, BJP, SUM, TMP;
	int J, JSUM, M;
	double BESSJ( int N, double X );
	double BESSJ0( double X );
	double BESSJ1( double X );
	TMP = 0.0;
	if( N == 0 ) return BESSJ0( X );
	if( N == 1 ) return BESSJ1( X );
	if( X == 0.0 ) return 0.0;
	TOX = 2.0 / X;
	if( X > N )
	{
		BJM = BESSJ0( X );
		BJ  = BESSJ1( X );
		for( J = 1; J < N; J++ )
		{
			BJP = J * TOX * BJ - BJM;
			BJM = BJ;
			BJ  = BJP;
			TMP = BJ;
		}
		return TMP;
	}
	else
	{
		M = ( int )( 2 * ( ( N + floor( sqrt( IACC * N ) ) ) / 2 ) );
		TMP = 0.0;
		JSUM = 0;
		SUM = 0.0;
		BJP = 0.0;
		BJ  = 1.0;
		for( J = M; J > 0; J-- )
		{
			BJM = J * TOX * BJ - BJP;
			BJP = BJ;
			BJ  = BJM;
			if( fabs( BJ ) > BIGNO )
			{
				BJ  *= BIGNI;
				BJP *= BIGNI;
				TMP *= BIGNI;
				SUM *= BIGNI;
			}
			if( JSUM != 0 )  SUM += BJ;
			JSUM = 1 - JSUM;
			if( J == N ) TMP = BJP;
			SUM = 2.0 * SUM - BJ;
		}
		return ( TMP / SUM );
	}
}
// ----------------------------------------------------------------------
double BESSJ0( double X )
{
	/*    THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
	      OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
	      BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
	      REFERENCES :
	      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	      VOL.5, 1962.  */
	double AX, FR, FS, Z, FP, FQ, TMP, XX;
	double Y, P1, P2, P3, P4, P5, R1, R2, R3, R4, R5, R6, Q1, Q2, Q3, Q4, Q5, S1, S2, S3, S4, S5, S6;
	P1 = 1.0; P2 = -0.1098628627E-2; P3 = 0.2734510407E-4;
	P4 = -0.2073370639E-5; P5 = 0.2093887211E-6;
	Q1 = -0.1562499995E-1; Q2 = 0.1430488765E-3; Q3 = -0.6911147651E-5;
	Q4 = 0.7621095161E-6; Q5 = -0.9349451520E-7;
	R1 = 57568490574.0; R2 = -13362590354.0; R3 = 651619640.7;
	R4 = -11214424.18; R5 = 77392.33017; R6 = -184.9052456;
	S1 = 57568490411.0; S2 = 1029532985.0; S3 = 9494680.718;
	S4 = 59272.64853; S5 = 267.8532712; S6 = 1.0;
	if( X == 0.0 ) return 1.0;
	AX = fabs( X );
	if( AX < 8 )
	{
		Y = X * X;
		FR = R1 + Y * ( R2 + Y * ( R3 + Y * ( R4 + Y * ( R5 + Y * R6 ) ) ) );
		FS = S1 + Y * ( S2 + Y * ( S3 + Y * ( S4 + Y * ( S5 + Y * S6 ) ) ) );
		TMP = FR / FS;
	}
	else
	{
		Z = 8. / AX;
		Y = Z * Z;
		XX = AX - 0.785398164;
		FP = P1 + Y * ( P2 + Y * ( P3 + Y * ( P4 + Y * P5 ) ) );
		FQ = Q1 + Y * ( Q2 + Y * ( Q3 + Y * ( Q4 + Y * Q5 ) ) );
		TMP = sqrt( 0.636619772 / AX ) * ( FP * cos( XX ) - Z * FQ * sin( XX ) );
	}
	return TMP;
}
// ----------------------------------------------------------------------
double BESSJ1( double X )
{
	/*    THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
	      OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
	      BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
	      REFERENCES :
	      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	      VOL.5, 1962.  */
	double AX, FR, FS, Z, FP, FQ, TMP, XX;
	double Y, P1, P2, P3, P4, P5, P6, R1, R2, R3, R4, R5, R6, Q1, Q2, Q3, Q4, Q5, S1, S2, S3, S4, S5, S6;
	P1 = 1.0; P2 = 0.183105E-2; P3 = -0.3516396496E-4;
	P4 = 0.2457520174E-5; P5 = -0.240337019E-6; P6 = 0.636619772;
	Q1 = 0.04687499995; Q2 = -0.2002690873E-3; Q3 = 0.8449199096E-5;
	Q4 = -0.88228987E-6; Q5 = 0.105787412E-6;
	R1 = 72362614232.0; R2 = -7895059235.0; R3 = 242396853.1;
	R4 = -2972611.439; R5 = 15704.48260; R6 = -30.16036606;
	S1 = 144725228442.0; S2 = 2300535178.0; S3 = 18583304.74;
	S4 = 99447.43394; S5 = 376.9991397; S6 = 1.0;
	AX = fabs( X );
	if( AX < 8 )
	{
		Y = X * X;
		FR = R1 + Y * ( R2 + Y * ( R3 + Y * ( R4 + Y * ( R5 + Y * R6 ) ) ) );
		FS = S1 + Y * ( S2 + Y * ( S3 + Y * ( S4 + Y * ( S5 + Y * S6 ) ) ) );
		TMP = X * ( FR / FS );
	}
	else
	{
		Z = 8. / AX;
		Y = Z * Z;
		XX = AX - 2.35619491;
		FP = P1 + Y * ( P2 + Y * ( P3 + Y * ( P4 + Y * P5 ) ) );
		FQ = Q1 + Y * ( Q2 + Y * ( Q3 + Y * ( Q4 + Y * Q5 ) ) );
		TMP = sqrt( P6 / AX ) * ( cos( XX ) * FP - Z * sin( XX ) * FQ ) * Sign( S6, X );
	}
	return TMP;
}

// end of file trootj.cpp

