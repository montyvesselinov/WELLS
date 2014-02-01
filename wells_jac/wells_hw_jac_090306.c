/* WELLSUB.C - Wells functions
 *             Velimir V. Vesselinov - Monty SP (c) 1992
 */
#include <math.h>
#include <stdio.h>

#include "design.h"
#include "wells2d2.h"

/* Prototypes of the functions in this file */
void  Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd , double *t_tot_obs, double *t_tot_prod );
void  Save_Problem( struct Problem_Info pi, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd );
void  Save_Well_S( char *fileName, struct Well_Data *wd );
void  Save_Point_S( char *fileName, struct Point_Data *pd );
void  Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd, struct Jacob_Info *ji );
void  Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd, struct Jacob_Info *ji );
/*void  Calc_Save_Grid( double time, struct Aquifer_Data ad, struct Well_Data *wd, struct Grid_Data *gd, char *fileName );*/
void  Jacob_Out( struct Well_Data *wd, struct Point_Data *pd, struct Jacob_Info *ji );
double Theis( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd , double *t_Jac , int * n_Jac );
double Theis_unc( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd, double *t_Jac, int *n_Jac);
double Hantush( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd, double *t_Jac, int *n_Jac );
double Hantush_unc( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd, double *t_Jac, int *n_Jac );
double Ei( double u );
double JacobsEi( double u );
double W( double ra, double rb );
double K0( double x );
int matherr( struct exception *except );

int main( int argn, char *args[] )
{
	struct Problem_Info pit;
	struct Aquifer_Data adt;
	struct Well_Data wdt;
	struct Point_Data pdt;
	struct Jacob_Info jit;

        if( argn < 4 )
        {
                printf( "Usage: wells in_problem_file out_well_drawdown out_point_drawdown\n" ); 
                exit( 1 );
        }

	strcpy( pit.file, args[1] );
	Load_Problem( &pit, &adt, &wdt, &pdt , &jit.t_tot_obs[0] , &jit.t_tot_prod[0] );
	strcpy( pit.file, "debug.wel" );
	Save_Problem( pit, adt, &wdt, &pdt );
	Calc_Wells( adt, &wdt , &jit );
	Save_Well_S( args[2], &wdt );
	Calc_Points( adt, &wdt, &pdt , &jit );
	Save_Point_S( args[3], &pdt );
	Jacob_Out( &wdt , &pdt , &jit );
}

void Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd, double *t_tot_obs, double *t_tot_prod )
{
	FILE *infile;
	double del_t;
	char buf[50];
	int  i, j;

	if ( ( infile = fopen( (*pi).file, "r" ) ) == NULL )
	{
		printf( "Problem file could not be read!\n" );
		return;
	}
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ":%[^\n]60s\n", (*pi).name );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &(*ad).co_aquifer );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i %i %f %f %f\n", &(*ad).co_x_axis, &(*ad).co_y_axis, &(*ad).origin_x, &(*ad).origin_y, &(*ad).origin_angle );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &(*wd).nW );
	for( i = 0; i < (*wd).nW; i++ )
	{
		fscanf( infile, "%s %lf %lf %lf %i\n", &(*wd).id[i], &(*wd).x[i], &(*wd).y[i], &(*wd).r[i], &(*wd).nQ[i] );
/*		printf( "%s %f %f %f %i\n", (*wd).id[i], (*wd).x[i], (*wd).y[i], (*wd).r[i], (*wd).nQ[i] ); */
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*wd).h[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*wd).m[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*wd).k[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*wd).S[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*wd).B[i] );
		(*wd).B[i] = pow( 10, (*wd).B[i] );
		for( j = 0; j < (*wd).nQ[i]; j++ )
			fscanf( infile, "%lf %lf\n", &(*wd).t[i][j], &(*wd).Q[i][j] );
		t_tot_prod[i] = (*wd).t[i][(*wd).nQ[i] - 1 ] - (*wd).t[i][0];
	}
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &(*pd).nP );
	for( i = 0; i < (*pd).nP; i++ )
	{
		fscanf( infile, "%s %lf %lf %i\n", &(*pd).id[i], &(*pd).x[i], &(*pd).y[i], &(*pd).nT[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &(*pd).h[i] );
		for( j = 0; j < (*pd).nT[i]; j++ )
			fscanf( infile, "%lf\n", &(*pd).t[i][j] );
		t_tot_obs[i] = (*pd).t[i][ (*pd).nT[i] - 1 ] - (*pd).t[i][0];
	}
	fclose( infile );
	for ( i = 0; i < (*pd).nP; i++ )
	{
		for ( j = 0; j < (*wd).nW; j++ )
		{
			(*ad).kobs[i][j] = pow( 10, (*wd).k[j] );
/*			printf( "%g %g %g\n", (*ad).kobs[i][j], (*wd).k[j], (*pd).k[i] ); */
			(*ad).Sobs[i][j] = pow( 10, (*wd).S[j] );
/*			printf( "%g %g %g\n", (*ad).Sobs[i][j], (*wd).S[j], (*pd).S[i] ); */
		}
	}
	for ( i = 0; i < (*wd).nW; i++ )
	{
		for ( j = 0; j < (*wd).nW; j++ )
		{
			(*ad).kprod[j][i] = pow( 10, ( (*wd).k[j] + (*wd).k[i] ) / 2 );
			(*ad).Sprod[j][i] = pow( 10, ( (*wd).S[j] + (*wd).S[i] ) / 2 );
		}
	}

}

void Save_Problem( struct Problem_Info pi, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd )
{
	FILE *outfile;
	int  i, j;

	if ( ( outfile = fopen( pi.file, "w" ) ) == NULL )
	{
		printf( "Problem file could not be opened!\n" );
		return;
	}
	fprintf( outfile, "Problem name: %s\n", pi.name );
	fprintf( outfile, "Aquifer type        : %i\n", ad.co_aquifer );
	fprintf( outfile, "Boundary information: %i %i %g %g %g\n", ad.co_x_axis, ad.co_y_axis, ad.origin_x, ad.origin_y, ad.origin_angle );
	fprintf( outfile, "--- Number of  wells: %i\n", (*wd).nW );
	for( i = 0; i < (*wd).nW; i++ )
	{
		fprintf( outfile, "%s %g %g %g %i\n", (*wd).id[i], (*wd).x[i], (*wd).y[i], (*wd).r[i], (*wd).nQ[i] );
		fprintf( outfile, "Initial Head        : %g\n", (*wd).h[i] );
		fprintf( outfile, "Aquifer Tickness    : %g\n", (*wd).m[i] );
		fprintf( outfile, "Permeability        : %g\n", (*wd).k[i] );
		fprintf( outfile, "Storage coefficient : %g\n", (*wd).S[i] );
		fprintf( outfile, "Leakage coefficient : %g\n", (*wd).B[i] );
		for( j = 0; j < (*wd).nQ[i]; j++ )
			fprintf( outfile, "%g %g\n", (*wd).t[i][j], (*wd).Q[i][j] );
	}
	fprintf( outfile, "--- Number of points: %i\n", (*pd).nP );
	for( i = 0; i < (*pd).nP; i++ )
	{
		fprintf( outfile, "%s %g %g %g %i\n", (*pd).id[i], (*pd).x[i], (*pd).y[i], (*pd).h[i], (*pd).nT[i] );
		fprintf( outfile, "Initial Head        : %g\n", (*pd).h[i] );
		fprintf( outfile, "Permeability        : %g\n", (*pd).k[i] );
		fprintf( outfile, "Storage coefficient : %g\n", (*wd).S[i] );
		for( j = 0; j < (*pd).nT[i]; j++ )
			fprintf( outfile, "%g\n", (*pd).t[i][j] );
	}
	fprintf(outfile, "Prod well, obs well, Harm. mean perm., geom. mean storage coeff.\n");
	for ( i = 0; i < (*pd).nP; i++ )
	{
		for ( j = 0; j < (*wd).nW; j++ )
		{
			fprintf(outfile, "%i %i %g %g\n", i, j, ad.kobs[i][j], ad.Sobs[i][j] );
		}
	}
	fprintf(outfile, "Prod well, prod well, Harm. mean perm., geom. mean storage coeff.\n");
	for ( i = 0; i < (*wd).nW; i++ )
	{
		for ( j = 0; j < (*wd).nW; j++ )
		{
			fprintf(outfile, "%i %i %g %g\n", i, j, ad.kprod[i][j], ad.Sprod[i][j] );
		}
	}

	fclose( outfile );
}

void Save_Well_S( char *fileName, struct Well_Data *wd )
{
	FILE *outfile;
	int i, j;

	if ( ( outfile = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Well drawdowns could not be saved!\n" );
		return;
	}
	for( i = 0; i < (*wd).nW; i++ )
	{
		fprintf( outfile, "%s %g %g %i\n", (*wd).id[i], (*wd).x[i], (*wd).y[i], (*wd).nQ[i] );
		for( j = 0; j < (*wd).nQ[i]; j++ )
			fprintf( outfile, "%.7g %g %g %.6f\n", (*wd).t[i][j], (*wd).Q[i][j], (*wd).s[i][j], (*wd).h[i] - (*wd).s[i][j] );
	}
	fclose( outfile );
}

void Save_Point_S( char *fileName, struct Point_Data *pd )
{
	FILE *outfile;
	int i, j;

	if ( ( outfile = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Data could not be saved!\n" );
		return;
	}
	for( i = 0; i < (*pd).nP; i++ )
	{
		fprintf( outfile, "%s %g %g %i\n", (*pd).id[i], (*pd).x[i], (*pd).y[i], (*pd).nT[i] );
		for( j = 0; j < (*pd).nT[i]; j++ )
			fprintf( outfile, "%.7g %g %.6f\n", (*pd).t[i][j], (*pd).s[i][j], (*pd).h[i] - (*pd).s[i][j] );
	}
	fclose( outfile );
}

void Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd , struct Jacob_Info *ji )
{
	int   i, j, code1, code2, code3;
	double x, y, xc, yc, ca, sa, h, k, S, del_t,
		  (*fp)( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd, double *t_Jac_prod, int *n_Jac_prod );

	switch( ad.co_aquifer )
	{
		case CONFINED:
			fp = Theis;
			break;
		case UNCONFINED:
			fp = Theis_unc;
			break;
		case LEAKY:
			fp = Hantush;
			break;
		case LEAKY_UNCONFINED:
			fp = Hantush_unc;
			break;
	}
	if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
	{
		ca = cos( ad.origin_angle );
		sa = sin( ad.origin_angle );
		code1 = code2 = code3 = 0;
		if( ad.co_x_axis != NO_BOUNDARY )
			code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_y_axis != NO_BOUNDARY )
			code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY )
			code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	}

	for( i = 0; i < (*wd).nW; i++ )
	{
		for( j = 0; j < (*wd).nW; j++ ) 
		{
			(*ji).n_Jac_prod[i][j] = 0;
			(*ji).t_Jac_prod[i][j] = 0;
		}
		for( j = 0; j < (*wd).nQ[i]; j++ )
		{
			if( j == 0 )
			{
				del_t = 0;
			} else {
				del_t = (*wd).t[i][j] - (*wd).t[i][j-1];
			}
			(*wd).s[i][j] = fp( (*wd).x[i], (*wd).y[i], (*wd).t[i][j], del_t, (*wd).h[i], ad.kprod[i], ad.Sprod[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
			if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
			{
				xc = (*wd).x[i] - ad.origin_x;
				yc = (*wd).y[i] - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
				if( code1 )
				{
					xc =x * ca - -y * sa + ad.origin_x;
					yc = -y * ca +x * sa + ad.origin_y;
					(*wd).s[i][j] += code1 * fp( xc, yc, (*wd).t[i][j], del_t, (*wd).h[i], ad.kprod[i], ad.Sprod[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
				if( code2 )
				{
					xc = -x * ca -y * sa + ad.origin_x;
					yc =y * ca + -x * sa + ad.origin_y;
					(*wd).s[i][j] += code2 * fp( xc, yc, (*wd).t[i][j], del_t, (*wd).h[i], ad.kprod[i], ad.Sprod[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					(*wd).s[i][j] += code3 * fp( xc, yc, (*wd).t[i][j], del_t, (*wd).h[i], ad.kprod[i], ad.Sprod[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
			}
		}
	}
}

void Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd , struct Jacob_Info *ji )
{
	int   i, j, code1, code2, code3;
	double x, y, xc, yc, ca, sa, h, k, S, del_t,
		  (*fp)( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd, double *t_Jac_prod, int *n_Jac_prod );

	switch( ad.co_aquifer )
	{
		case CONFINED:
			fp = Theis;
			break;
		case UNCONFINED:
			fp = Theis_unc;
			break;
		case LEAKY:
			fp = Hantush;
			break;
		case LEAKY_UNCONFINED:
			fp = Hantush_unc;
			break;
	}
	if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
	{
		ca = cos( ad.origin_angle );
		sa = sin( ad.origin_angle );
		code1 = code2 = code3 = 0;
		if( ad.co_x_axis != NO_BOUNDARY )
			code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_y_axis != NO_BOUNDARY )
			code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
		if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY )
			code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	}

	for( i = 0; i < (*pd).nP; i++ )
	{
		for( j = 0; j < (*wd).nW; j++ ) 
		{
			(*ji).n_Jac_obs[i][j] = 0;
			(*ji).t_Jac_obs[i][j] = 0;
		}
		for( j = 0; j < (*pd).nT[i]; j++ )
		{
			if( j == 0 )
			{
				del_t = 0;
			} else {
				del_t = (*pd).t[i][j] - (*pd).t[i][j-1];
			}
			(*pd).s[i][j] = fp( (*pd).x[i], (*pd).y[i], (*pd).t[i][j], del_t , (*pd).h[i], ad.kobs[i], ad.Sobs[i], wd , &(*ji).t_Jac_obs[i][0] , &(*ji).n_Jac_obs[i][0] );
			if( ad.co_x_axis != NO_BOUNDARY || ad.co_y_axis != NO_BOUNDARY )
			{
				xc = (*pd).x[i] - ad.origin_x;
				yc = (*pd).y[i] - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
				if( code1 )
				{
					xc =x * ca - -y * sa + ad.origin_x;
					yc = -y * ca +x * sa + ad.origin_y;
					(*pd).s[i][j] += code1 * fp( xc, yc, (*pd).t[i][j], del_t, (*pd).h[i], ad.kobs[i], ad.Sobs[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
				if( code2 )
				{
					xc = -x * ca -y * sa + ad.origin_x;
					yc =y * ca + -x * sa + ad.origin_y;
					(*pd).s[i][j] += code2 * fp( xc, yc, (*pd).t[i][j], del_t, (*pd).h[i], ad.kobs[i], ad.Sobs[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					(*pd).s[i][j] += code3 * fp( xc, yc, (*pd).t[i][j], del_t, (*pd).h[i], ad.kobs[i], ad.Sobs[i], wd , &(*ji).t_Jac_prod[i][0] , &(*ji).n_Jac_prod[i][0] );
				}
			}
		}
	}
}
/* Warning: ad.h no longer valid
void Calc_Save_Grid( double time, struct Aquifer_Data ad, struct Well_Data *wd, struct Grid_Data *gd, char *fileName )
{
	FILE  *out;
	int   i, j, code1, code2, code3;
	char  buf[30];
	double s, xg, yg, dx, dy, x, y, xc, yc, ca, sa, h,
		  (*fp)( double x, double y, double time, double h, struct Well_Data *wd );

	if ( ( out = fopen( fileName, "w" ) ) == NULL )
	{
		printf( "Grid file could not be opened!" );
		return;
	}
	fprintf( out, "DSAA\n" );
	fprintf( out, "%i %i\n", (*gd).nx, (*gd).ny );
	fprintf( out, "%.3g %.3g\n", (*gd).xmin, (*gd).xmax );
	fprintf( out, "%.3g %.3g\n", (*gd).ymin, (*gd).ymax );
	fprintf( out, "%.3g %.3g\n", 0.0, 0.0 );
	switch( ad.co_aquifer )
	{
		case CONFINED:
			fp = Theis;
			break;
		case UNCONFINED:
			fp = Theis_unc;
			break;
		case LEAKY:
			fp = Hantush;
			break;
		case LEAKY_UNCONFINED:
			fp = Hantush_unc;
			break;
	}

	dx = ( (*gd).xmax - (*gd).xmin ) / ( (*gd).nx - 1 );
	dy = ( (*gd).ymax - (*gd).ymin ) / ( (*gd).ny - 1 );
	ca = cos( ad.origin_angle );
	sa = sin( ad.origin_angle );
	code1 = code2 = code3 = 0;
	if( ad.co_x_axis != NO_BOUNDARY )
		code1 = ( ad.co_x_axis == TYPE_1 ) ? -1 : 1;
	if( ad.co_y_axis != NO_BOUNDARY )
		code2 = ( ad.co_y_axis == TYPE_1 ) ? -1 : 1;
	if( ad.co_x_axis != NO_BOUNDARY && ad.co_y_axis != NO_BOUNDARY )
		code3 = ( ad.co_x_axis != ad.co_y_axis ) ? -1 : 1;
	for( yg = (*gd).ymin, j = 0; j < (*gd).ny; j++, yg += dy )
	{
		printf( "\rRow %4i of %4i", j + 1, (*gd).ny );
		for( xg = (*gd).xmin, i = 0; i < (*gd).nx; i++, xg += dx, fprintf( out, "%12g ", s ) )
		{
			s = fp( xg, yg, time, ad.h, wd );
			if( code1 || code2 )
			{
				xc = xg - ad.origin_x;
				yc = yg - ad.origin_y;
				x = xc * ca + xc * sa;
				y = yc * ca - yc * sa;
			}
			else
				continue;
			if( code1 )
			{
				xc =  x * ca - -y * sa + ad.origin_x;
				yc = -y * ca +  x * sa + ad.origin_y;
				s += code1 * fp( xc, yc, time, ad.h, wd );
			}
			if( code2 )
			{
				xc = -x * ca -  y * sa + ad.origin_x;
				yc =  y * ca + -x * sa + ad.origin_y;
				s += code2 * fp( xc, yc, time, ad.h, wd );
			}
			if( code3 )
			{
				xc = -x * ca - -y * sa + ad.origin_x;
				yc = -y * ca + -x * sa + ad.origin_y;
				s += code3 * fp( xc, yc, time, ad.h, wd );
			}
		}
		fprintf( out, "\n" );
	}
	fclose( out );
}
*/
double Theis( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd , double *t_Jac , int * n_Jac )
{
	double dx, dy, r, T, sub1, sub2, sub3, sub4, s;
	int i, n, nR;

	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		T = (*wd).m[i] * k[i]; /* transmissivity */
/*		printf( "%i %g %g\n", i, k[i], S[i] ); */
		sub1 = 12.5663706 * T; /* 4 * PI = 12.5663706 */
		sub2 = S[i] / ( 4.0 * T );
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = dx * dx + dy * dy; /* squared radial distance */
		sub3 = (*wd).r[i] * (*wd).r[i]; /* squared well radius */
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		sub4 = time - (*wd).t[i][0];
		if( sub4 > 0 )
			s += (*wd).Q[i][0] * JacobsEi( sub3 / sub4 ) / sub1; /* initial drawdown */
		nR = (*wd).nQ[i]; /* number of step pumping rate changes */
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ ) /* loop through all the pumping rates */
		{
			sub4 = time - (*wd).t[i][n];
			if( sub4 > 0 )
				s += ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * JacobsEi( sub3 / sub4 ) / sub1;
		}
		if( sub3 / sub4 < 0.03 ) 
		{
			t_Jac[i] += del_t;
			n_Jac[i]++;
		}
	}
	return( s );
}

double Theis_unc( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd , double *t_Jac , int * n_Jac )
{
	double dx, dy, r, h2, sub1, sub2, sub3, sub4, s;
	int i, n, nR;

	h2 = h * h;
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		sub1 = 6.2831853 * k[i];
		sub2 = S[i] / ( 4.0 * k[i] * (*wd).m[i] );
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = dx * dx + dy * dy;
		sub3 = (*wd).r[i] * (*wd).r[i];
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		sub4 = time - (*wd).t[i][0];
		if( sub4 > 0 )
		{
			s += h - sqrt( h2 - (*wd).Q[i][0] * Ei( sub3 / sub4 ) / sub1 );
		}
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
		{
			sub4 = time - (*wd).t[i][n];
			if( sub4 > 0 )
				s += h - sqrt( h2 - ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * JacobsEi( sub3 / sub4 ) / sub1 );
		}
		if( sub3 / sub4 < 0.03 ) 
		{
			t_Jac[i] += del_t;
			n_Jac[i]++;
		}

	}
	return( s );
}

double Hantush( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd , double *t_Jac , int * n_Jac )
{
	double dx, dy, r, T, sub1, sub2, sub3, sub4, sub5, s;
	int i, n, nR;

	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		T = (*wd).m[i] * k[i];
		sub1 = 12.5663706 * T;
		sub2 = S[i] / ( 4.0 * T );
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < (*wd).r[i] ) r = (*wd).r[i];
		sub3 = r * r * sub2;
		sub4 = r / (*wd).B[i];
		sub5 = time - (*wd).t[i][0];
		if( sub5 > 0 )
			s += (*wd).Q[i][0] * W( sub3 / sub5, sub4 ) / sub1;
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
		{
			sub5 = time - (*wd).t[i][n];
			if( sub5 > 0 )
				s += ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * W( sub3 / sub5, sub4 ) / sub1;
		}
	}
	return( s );
}

double Hantush_unc( double x, double y, double time, double del_t, double h, double *k, double *S, struct Well_Data *wd , double *t_Jac , int * n_Jac )
{
	double dx, dy, r, h2, sub1, sub2, sub3, sub4, sub5, s;
	int i, n, nR;

	h2 = h * h;
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		sub1 = 6.2831853 * k[i];
		sub2 = S[i] / ( 4.0 * k[i] * (*wd).m[i] );
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < (*wd).r[i] ) r = (*wd).r[i];
		sub3 = r * r * sub2;
		sub4 = r / (*wd).B[i];
		sub5 = time - (*wd).t[i][0];
		if( sub5 > 0 )
			s += h - sqrt( h2 - (*wd).Q[i][0] * W( ( sub3 / sub5 ), sub4 ) / sub1 );
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
		{
			sub5 = time - (*wd).t[i][n];
			if( sub5 > 0 )
				s += h - sqrt( h2 - ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * W( ( sub3 / sub5 ), sub4 ) / sub1 );
		}
	}
	return( s );
}

/*  The Theis Function - Ei(x)
 */
double Ei( double u )
{
	double s1, s2;

	if ( u > 60.0 || u <= 0.0 ) return( 0.0 );
	if ( u < 1.0 )
	{
		s1 = ( ( ( ( 0.00107857 * u - 0.00976004 ) * u + 0.05519968 ) * u - 0.24991055 ) * u + 0.99999193 ) * u - 0.57721566;
		return( s1 - log ( u ) );
	}
	else
	{
		s1 = ( ( ( u + 8.57332875 ) * u + 18.05901697 ) * u + 8.63476088 ) * u + 0.267773734;
		s2 = ( ( ( u + 9.57332875 ) * u + 25.63295615 ) * u + 21.09965308 ) * u + 3.958496923;
		return( ( 1 / u ) * ( s1 / s2 ) * exp( -u ) );
	}
}
double JacobsEi( double u )
{
	return( -0.57721566 - log( u ) );
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

void  Jacob_Out( struct Well_Data *wd, struct Point_Data *pd, struct Jacob_Info *ji )
{
	int i, j;
	double perc_t, perc_n;
	printf("Percentage of calculations and time with valid Jacob's approximation (u < 0.03)\n");
	printf("Point Well %_of_Calcs %_of_Time\n");
	for( i = 0; i < (*pd).nP; i++ )
	{
		for( j = 0; j < (*wd).nW; j++ )
		{
			perc_t = (*ji).t_Jac_obs[i][j] / (*ji).t_tot_obs[i];
			perc_n = (double) (*ji).n_Jac_obs[i][j] / (*pd).nT[i];
			printf("%s %s %f %f\n", (*pd).id[i], (*wd).id[j], perc_n, perc_t );
			//printf("%s %s %d %d %6.1f %6.1f\n", (*pd).id[i], (*wd).id[j], (*pd).nT[i], (*ji).n_Jac_obs[i][j], (*ji).t_tot_obs[i], (*ji).t_Jac_obs[i][j]);
		}
	}
	for( i = 0; i < (*wd).nW; i++ )
	{
		for( j = 0; j < (*wd).nW; j++ )
		{
			perc_t = (*ji).t_Jac_prod[i][j] / (*ji).t_tot_prod[i];
			perc_n = (double) (*ji).n_Jac_prod[i][j] / (*wd).nQ[i];
			printf("%s %s %f %f\n", (*wd).id[i], (*wd).id[j], perc_n, perc_t );
			//printf("%s %s %d %d %6.1f %6.1f\n", (*wd).id[i], (*wd).id[j], (*wd).nQ[i], (*ji).n_Jac_prod[i][j], (*ji).t_tot_prod[i], (*ji).t_Jac_prod[i][j]);
		}
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
	}
	else
	{
		printf( "ERROR!" );
		return 0;
	}
}
