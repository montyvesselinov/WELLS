/* WELLSUB.C - Wells functions
 *             Velimir V. Vesselinov - Monty SP (c) 1992
 */
#include <math.h>
#include <stdio.h>

#include "design.h"
#include "wells.h"

/* Prototypes of the functions in this file */
void  Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd );
void  Save_Problem( struct Problem_Info pi, struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd );
void  Save_Well_S( char *fileName, struct Aquifer_Data ad, struct Well_Data *wd );
void  Save_Point_S( char *fileName, struct Aquifer_Data ad, struct Point_Data *pd );
void  Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd );
void  Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd );
void  Calc_Save_Grid( float time, struct Aquifer_Data ad, struct Well_Data *wd, struct Grid_Data *gd, char *fileName );
float Theis( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );
float Theis_unc( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );
float Hantush( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );
float Hantush_unc( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );
float Ei( float u );
float W( float ra, float rb );
float K0( float x );
int matherr( struct exception *except );

int main( int argn, char *args[] )
{
	struct Problem_Info pit;
	struct Aquifer_Data adt;
	struct Well_Data wdt;
	struct Point_Data pdt;

        if( argn < 4 )
        {
                printf( "Usage: wells in_problem_file out_well_drawdown out_point_drawdown\n" ); 
                exit( 1 );
        }

	strcpy( pit.file, args[1] );
	Load_Problem( &pit, &adt, &wdt, &pdt );
	strcpy( pit.file, "debug.wel" );
	Save_Problem( pit, adt, &wdt, &pdt );
	Calc_Wells( adt, &wdt );
	Save_Well_S( args[2], adt, &wdt );
	Calc_Points( adt, &wdt, &pdt );
	Save_Point_S( args[3], adt, &pdt );
}

void Load_Problem( struct Problem_Info *pi, struct Aquifer_Data *ad, struct Well_Data *wd, struct Point_Data *pd )
{
	FILE *infile;
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
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %f\n", &(*ad).h );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %f\n", &(*ad).m );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %f\n", &(*ad).k );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %f\n", &(*ad).S );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %f\n", &(*ad).B );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &(*wd).nW );
	for( i = 0; i < (*wd).nW; i++ )
	{
		fscanf( infile, "%s %f %f %f %i\n", &(*wd).id[i], &(*wd).x[i], &(*wd).y[i], &(*wd).r[i], &(*wd).nQ[i] );
		for( j = 0; j < (*wd).nQ[i]; j++ )
			fscanf( infile, "%f %f\n", &(*wd).t[i][j], &(*wd).Q[i][j] );
	}
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &(*pd).nP );
	for( i = 0; i < (*pd).nP; i++ )
	{
		fscanf( infile, "%s %f %f %i\n", &(*pd).id[i], &(*pd).x[i], &(*pd).y[i], &(*pd).nT );
		for( j = 0; j < (*pd).nT[i]; j++ )
			fscanf( infile, "%f\n", &(*pd).t[i][j] );
	}
	fclose( infile );
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
	fprintf( outfile, "Initial Head        : %g\n", ad.h );
	fprintf( outfile, "Aquifer Tickness    : %g\n", ad.m );
	fprintf( outfile, "Permeability        : %g\n", ad.k );
	fprintf( outfile, "Storage coefficient : %g\n", ad.S );
	fprintf( outfile, "Leakage coefficient : %g\n", ad.B );
	fprintf( outfile, "--- Number of  wells: %i\n", (*wd).nW );
	for( i = 0; i < (*wd).nW; i++ )
	{
		fprintf( outfile, "%s %g %g %g %i\n", (*wd).id[i], (*wd).x[i], (*wd).y[i], (*wd).r[i], (*wd).nQ[i] );
		for( j = 0; j < (*wd).nQ[i]; j++ )
			fprintf( outfile, "%g %g\n", (*wd).t[i][j], (*wd).Q[i][j] );
	}
	fprintf( outfile, "--- Number of points: %i\n", (*pd).nP );
	for( i = 0; i < (*pd).nP; i++ )
	{
		fprintf( outfile, "%s %g %g %i\n", (*pd).id[i], (*pd).x[i], (*pd).y[i], (*pd).nT[i] );
		for( j = 0; j < (*pd).nT[i]; j++ )
			fprintf( outfile, "%g\n", (*pd).t[i][j] );
	}
	fclose( outfile );
}

void Save_Well_S( char *fileName, struct Aquifer_Data ad, struct Well_Data *wd )
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
			fprintf( outfile, "%g %g %g %g\n", (*wd).t[i][j], (*wd).Q[i][j], (*wd).s[i][j], ad.h - (*wd).s[i][j] );
	}
	fclose( outfile );
}

void Save_Point_S( char *fileName, struct Aquifer_Data ad, struct Point_Data *pd )
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
			fprintf( outfile, "%g %g %g\n", (*pd).t[i][j], (*pd).s[i][j], ad.h - (*pd).s[i][j] );
	}
	fclose( outfile );
}

void Calc_Wells( struct Aquifer_Data ad, struct Well_Data *wd )
{
	int   i, j, code1, code2, code3;
	float x, y, xc, yc, ca, sa,
		  (*fp)( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );

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
		for( j = 0; j < (*wd).nQ[i]; j++ )
		{
			(*wd).s[i][j] = fp( (*wd).x[i], (*wd).y[i], (*wd).t[i][j], ad, wd );
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
					(*wd).s[i][j] += code1 * fp( xc, yc, (*wd).t[i][j], ad, wd );
				}
				if( code2 )
				{
					xc = -x * ca -y * sa + ad.origin_x;
					yc =y * ca + -x * sa + ad.origin_y;
					(*wd).s[i][j] += code2 * fp( xc, yc, (*wd).t[i][j], ad, wd );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					(*wd).s[i][j] += code3 * fp( xc, yc, (*wd).t[i][j], ad, wd );
				}
			}
		}
	}
}

void Calc_Points( struct Aquifer_Data ad, struct Well_Data *wd, struct Point_Data *pd )
{
	int   i, j, code1, code2, code3;
	float x, y, xc, yc, ca, sa,
		  (*fp)( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );

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
		for( j = 0; j < (*pd).nT[i]; j++ )
		{
			(*pd).s[i][j] = fp( (*pd).x[i], (*pd).y[i], (*pd).t[i][j], ad, wd );
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
					(*pd).s[i][j] += code1 * fp( xc, yc, (*pd).t[i][j], ad, wd );
				}
				if( code2 )
				{
					xc = -x * ca -y * sa + ad.origin_x;
					yc =y * ca + -x * sa + ad.origin_y;
					(*pd).s[i][j] += code2 * fp( xc, yc, (*pd).t[i][j], ad, wd );
				}
				if( code3 )
				{
					xc = -x * ca - -y * sa + ad.origin_x;
					yc = -y * ca + -x * sa + ad.origin_y;
					(*pd).s[i][j] += code3 * fp( xc, yc, (*pd).t[i][j], ad, wd );
				}
			}
		}
	}
}

void Calc_Save_Grid( float time, struct Aquifer_Data ad, struct Well_Data *wd, struct Grid_Data *gd, char *fileName )
{
	FILE  *out;
	int   i, j, code1, code2, code3;
	char  buf[30];
	float s, xg, yg, dx, dy, x, y, xc, yc, ca, sa,
		  (*fp)( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd );

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
			s = fp( xg, yg, time, ad, wd );
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
				s += code1 * fp( xc, yc, time, ad, wd );
			}
			if( code2 )
			{
				xc = -x * ca -  y * sa + ad.origin_x;
				yc =  y * ca + -x * sa + ad.origin_y;
				s += code2 * fp( xc, yc, time, ad, wd );
			}
			if( code3 )
			{
				xc = -x * ca - -y * sa + ad.origin_x;
				yc = -y * ca + -x * sa + ad.origin_y;
				s += code3 * fp( xc, yc, time, ad, wd );
			}
		}
		fprintf( out, "\n" );
	}
	fclose( out );
}

float Theis( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd )
{
	float dx, dy, r, T, sub1, sub2, sub3, s;
	int i, n, nR;

	T = ad.m * ad.k; /* transmissivity */
	sub1 = 12.5663706 * T; /* 4 * PI = 12.5663706 */
	sub2 = ad.S / ( 4.0 * T );
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = dx * dx + dy * dy; /* squared radial distance */
		sub3 = (*wd).r[i] * (*wd).r[i]; /* squared well radius */
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		s += (*wd).Q[i][0] * Ei( sub3 / ( time - (*wd).t[i][0] ) ) / sub1; /* initial drawdown */
		nR = (*wd).nQ[i]; /* number of step pumping rate changes */
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ ) /* loop through all the pumping rates */
			s += ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * Ei( sub3 / ( time - (*wd).t[i][n] ) ) / sub1;
	}
	return( s );
}

float Theis_unc( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd )
{
	float dx, dy, r, h2, sub1, sub2, sub3, s;
	int i, n, nR;

	h2 = ad.h * ad.h;
	sub1 = 6.2831853 * ad.k;
	sub2 = ad.S / ( 4.0 * ad.k * ad.m );
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = dx * dx + dy * dy;
		sub3 = (*wd).r[i] * (*wd).r[i];
		if( r < sub3 ) r = sub3;
		sub3 = r * sub2;
		s += ad.h - sqrt( h2 - (*wd).Q[i][0] * Ei( sub3 / ( time - (*wd).t[i][0] ) ) / sub1 );
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
			s += ad.h - sqrt( h2 - ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * Ei( sub3 / ( time - (*wd).t[i][n] ) ) / sub1 );
	}
	return( s );
}

float Hantush( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd )
{
	float dx, dy, r, T, sub1, sub2, sub3, sub4, s;
	int i, n, nR;

	T = ad.m * ad.k;
	sub1 = 12.5663706 * T;
	sub2 = ad.S / ( 4.0 * T );
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < (*wd).r[i] ) r = (*wd).r[i];
		sub3 = r * r * sub2;
		sub4 = r / ad.B;
		s += (*wd).Q[i][0] * W( sub3 / ( time - (*wd).t[i][0] ), sub4 ) / sub1;
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
			s += ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * W( sub3 / ( time - (*wd).t[i][n] ), sub4 ) / sub1;
	}
	return( s );
}

float Hantush_unc( float x, float y, float time, struct Aquifer_Data ad, struct Well_Data *wd )
{
	float dx, dy, r, h2, sub1, sub2, sub3, sub4, s;
	int i, n, nR;

	h2 = ad.h * ad.h;
	sub1 = 6.2831853 * ad.k;
	sub2 = ad.S / ( 4.0 * ad.k * ad.m );
	s = 0;
	for( i = 0; i < (*wd).nW; i++ )
	{
		dx = x - (*wd).x[i];	dy = y - (*wd).y[i];	r = sqrt( dx * dx + dy * dy );
		if( r < (*wd).r[i] ) r = (*wd).r[i];
		sub3 = r * r * sub2;
		sub4 = r / ad.B;
		s += ad.h - sqrt( h2 - (*wd).Q[i][0] * W( sub3 / ( time - (*wd).t[i][0] ), sub4 ) / sub1 );
		nR = (*wd).nQ[i];
		for( n = 1; n < nR && (*wd).t[i][n] < time; n++ )
			s += ad.h - sqrt( h2 - ( (*wd).Q[i][n] - (*wd).Q[i][n - 1] ) * W( sub3 / ( time - (*wd).t[i][n] ), sub4 ) / sub1 );
	}
	return( s );
}

/*  The Theis Function - Ei(x)
 */
float Ei( float u )
{
	float s1, s2;

	if ( u > 60.0 || u == 0.0 ) return( 0.0 );
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

/*  The Hantush Function
 *  The function is calculated by logarithmic integration using Simpson's rule
 *  source: W.Kinzelbach - Groundwater modelling
 */
float W( float ra, float rb )
{
	float ug, og, hi, sub1, sub2, x2, x4, s2, s4;
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
float K0( float x )
{
	float x1, x2;
	float s1, s2;

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
	}
	else
	{
		printf( "ERROR!" );
		return 0;
	}
}
