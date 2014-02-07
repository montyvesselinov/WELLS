#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define MAX_REC_SIZE 80
#define MAXWELL 20
#define MAXRATE 50000
#define MAXDELT 1000000

double CJ_Cutoff_Brents(double max_del_t, double T1, double T2, double S1, double S2, double r);
double CJ_Cutoff_Newton( double init_val, double T1, double T2, double S1, double S2, double r);
double Well_Func_u ( double t, void *params);
double Well_Func_u_deriv ( double t, void *params );
void Well_Func_u_fdf ( double t, void *params, double *y, double *dy);
struct u_params { double T1; double T2; double S1; double S2; double r; };
 
main( int argn, char *args[] )
{
	FILE *infile,*outfile,*tmpfile;
	int i,j,k,l,nW,nQ[20],nbins[20],bin[100],cal_start,Qp;
	char file[50],*bsnm,plotfl[50],buf[50],probname[50],id[20],pnt_id[20];
	double xp,yp,x,y,r,t[MAXRATE],del_t,max_del_t,cutoff[100],bin_avg,dist,cal_pt[MAXRATE];
	double T1,T2,S1,S2,cj_exp_t_brent, cj_exp_t_newton, cj_t, cj_init;

	strcpy(file, args[1] );

//	strcpy( buf, args[2] );
//	sscanf( buf, "%lf", &cj_init );
//	printf("Initial guess for Cooper-Jacob cutoff time: %lf\n", cj_init);
	
	for( i = 2; i < argn; i++ )
	{
		strcpy(buf, args[i] );
		sscanf(buf, "%i", &nbins[i-2] );
	}

	strcpy(buf, "awk '{if( $4 ~ /points/ ) print NR}' ");
	strcat(buf, file);
	strcat(buf, " > temp" );
	system(buf);

	if ( ( tmpfile = fopen( "temp", "r" ) ) == NULL )
	{
		printf( "Calibration line file could not be opened!\n" );
		return;
	}
	fscanf( tmpfile, "%i\n", &cal_start );	
	remove("temp");

	if ( ( infile = fopen( file, "r" ) ) == NULL )
	{
		printf( "WEL file could not be opened!\n" );
		return;
	}

	for( i = 1; i < cal_start + 1; i++ )
	{
		fgets( buf, MAX_REC_SIZE, infile );
	}
	fgets( buf, MAX_REC_SIZE, infile );
	sscanf( buf, "%s %lf %lf %i\n", &id, &x, &y, &Qp );	
	printf( "Number of observation times: %i\n", Qp );
	fgets( buf, MAX_REC_SIZE, infile );

	i = 0;
	while( fgets( buf, MAX_REC_SIZE, infile ) )
	{
		sscanf( buf, "%lf\n", &cal_pt[i] );
		i++;
	}
	fclose(infile);

	if ( ( infile = fopen( file, "r" ) ) == NULL )
	{
		printf( "WEL file could not be opened!\n" );
		return;
	}

	strcpy(buf, "sed -n '/points:/{n;p}' ");
	strcat(buf, file);
	strcat(buf, " > temp");
	system(buf); 

	if ( ( tmpfile = fopen( "temp", "r" ) ) == NULL )
	{
		printf( "Point coord file could not be opened!\n" );
		return;
	}

	fscanf( tmpfile, "%s %lf %lf %lf %i\n", &pnt_id, &xp, &yp, &r, &i );	
	printf( "Observation well name: %s\n", pnt_id );
	printf( "Observation x and y coords: %lf %lf\n", xp, yp );
	fclose(tmpfile);	
	remove("temp");

	for( i = 0; i < 3; i++ )
		fgets( buf, MAX_REC_SIZE, infile );

	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &nW );
	printf( "Number of wells: %d\n\n", nW);
	
	if( nW > argn - 2 )
	{
		printf( "Not enough arguments for %i wells!\n", nW );
		printf( "Usage: par-hist WELLFLNM #BINS1 #BINS2 ...\n");
		return;
	}	

	bsnm = strtok( file, ".");
	for( i = 0; i < nW; i++ )
	{
		fscanf( infile, "%s %lf %lf %lf %i\n", &id, &x, &y, &r, &nQ[i] );	
		printf( "Calculating histogram and Cooper-Jacob cutoff for %s/%s...\n", id, pnt_id);
		printf( "Number of rate changes: %i\n", nQ[i]);
		for( j = 0; j < 2; j++ )
			fgets( buf, MAX_REC_SIZE, infile );

		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf\n", &T1, &T2 );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf\n", &S1, &S2 );
		printf("%s/%s T = 10^%lf*exp(%lf/del_t)\n", id, pnt_id, T1, T2);
		printf("%s/%s S = 10^%lf*exp(%lf/del_t)\n", id, pnt_id, S1, S2);
		
		fgets( buf, MAX_REC_SIZE, infile );
		
		printf("Number of bins: %i\n", nbins[i]);
		
		for( j = 0; j < nQ[i]; j++ )
			fscanf( infile, "%lf %s\n", &t[j], buf );

		max_del_t = t[nQ[i]-1] - t[0];	
		printf("Max delta t: %lf days\n", max_del_t);

		for( j = 0; j < nbins[i] + 1; j++ )
		{
			cutoff[j] = max_del_t/nbins[i]*(j);
		}
		
		for(j = 0; j < nbins[i]; j++)
			bin[j] = 0;
	
		for(j = 0; j < Qp; j++ )
		{
			for(k = 0; k < nQ[i]; k++ )
			{
				if( cal_pt[j] > t[k] )
				{
					del_t = cal_pt[j] - t[k];
					for( l = 0; l < nbins[i]; l++ )
						if( del_t > cutoff[l] && del_t <= cutoff[l+1] )
						{
							bin[l]++;
							break;
						}
				} 
				else 
				{
					break;
				}
			}
		}

		dist = sqrt( pow(x - xp,2) + pow(y - yp,2));
		printf( "Distance between %s and %s: %lf\n", id, pnt_id, dist);

		strcpy(plotfl,bsnm);
		if( nW > 1 )
		{
			strcat( plotfl, "-" );
			strcat( plotfl, id);
			strcat( plotfl, ".gplot");
		} else {
			strcat( plotfl, ".gplot");
		}

		if ( ( outfile = fopen( plotfl, "w" ) ) == NULL )
		{
			printf( "Gnuplot temp file could not be opened!\n" );
			return;
		}
		
		if ( nbins[i] <= 11 )
			printf( "bin lb(d) ub(d) mid(d) count\n" );
		for(j = 0; j < nbins[i]; j++)
		{
			bin_avg = ( cutoff[j] + cutoff[j+1] ) / 2;
			if ( nbins[i] <= 11 )
				printf("%i %lf %lf %lf %i\n", j + 1, cutoff[j], cutoff[j + 1], bin_avg, bin[j]); 
			fprintf( outfile, "%i %lf %lf %lf %i\n", j + 1, cutoff[j], cutoff[j + 1], bin_avg, bin[j]); 
		}
		fclose(outfile);		
		printf("\nHistogram data stored in: %s\n", plotfl);
		printf( "\n" );
	
		printf( "Calculating Cooper-Jacob cutoff time...\n");	
		cj_t = pow(dist,2) * pow(10,S1) / (0.12 * pow(10,T1) );	
		cj_exp_t_brent = CJ_Cutoff_Brents( max_del_t, T1, T2, S1, S2, dist);
		cj_exp_t_newton = CJ_Cutoff_Newton( cj_t, T1, T2, S1, S2, dist);
		printf( "\nConstant Cooper-Jacob cutoff: %lf\n\n", cj_t);		

		FILE *pipe = popen("gnuplot -persist","w");
		fprintf(pipe, "set xlabel 'Time, days'\n");
		fprintf(pipe, "set ylabel 'Number of time interals'\n");
		strcpy(buf, "plot '");
		strcat(buf, plotfl);
		strcat(buf, "' us 4:5 with boxes\n");
		fprintf(pipe, buf);
		close(pipe);

	}	


	fclose( infile);

	return;

}

double CJ_Cutoff_Brents(double max_del_t, double T1, double T2, double S1, double S2, double r)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double rt=10;
	double x_lo = 1, x_hi = 10000, test_lo, test_hi;
	gsl_function F;
	struct u_params params = {T1, T2, S1, S2, r};	

	test_lo = Well_Func_u( x_lo, &params );
	test_hi = Well_Func_u( x_hi, &params );

	if ( test_lo * test_hi > 0 )
	{
		printf( "Brent's method stopped, y = 0 is not included in interval: %lf %lf!\n", test_lo, test_hi);
		return;
	}

	F.function = &Well_Func_u;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	printf( "Using %s method to identify exponential Cooper-Jacob cutoff...\n", gsl_root_fsolver_name (s));
	printf ("%5s [%9s, %9s] %9s\n", "iter", "lower", "upper", "root");

       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           rt = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi,
                                            0.001, .001);
     
           if (status == GSL_SUCCESS)
             printf ("Brent's method exponential Cooper-Jacob cutoff:\n");
     
           printf ("%5d [%f, %f] %f\n", iter, x_lo, x_hi, rt);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);
     
       return rt;

}

double CJ_Cutoff_Newton( double init_val, double T1, double T2, double S1, double S2, double r)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = init_val, test_deriv;
	gsl_function_fdf FDF;
	struct u_params params = { T1, T2, S1, S2, r};

	test_deriv = Well_Func_u_deriv( x, &params );

	if ( test_deriv == 0 )
	{
		printf( "Newtons method stopped, derivative is zero!\n\n");
		return;
	}

	FDF.f = &Well_Func_u;
	FDF.df = &Well_Func_u_deriv;
	FDF.fdf = &Well_Func_u_fdf;
	FDF.params = &params;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, x);

	printf( "\nUsing %s method to identify exponential Cooper-Jacob cutoff...\n", gsl_root_fdfsolver_name (s) );
	printf ("%-5s %10s\n", "iter", "root");

	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		test_deriv = Well_Func_u_deriv( x, &params );
		if ( test_deriv == 0 )
		{
			printf( "Newtons method stopped, derivative is zero!\n\n");
			return;
		}
		status = gsl_root_test_delta (x, x0, 1e-6, 1e-6);

		if (status == GSL_SUCCESS)
			printf( "Newton's method exponential Cooper-Jacob cutoff:\n" );

		printf( "%5d %10.7f\n", iter, x);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fdfsolver_free (s);
	
	return x;

}


double Well_Func_u ( double t, void *params)
{
	struct u_params *p = (struct u_params *) params;
	double ans;

	double T1 = p->T1;
	double T2 = p->T2;
	double S1 = p->S1;
	double S2 = p->S2;
	double r = p->r;

	ans = ( 0.03 - ( pow(r,2) * pow(10,S1) * exp(S2 / t) )/( 4 * pow(10,T1) * exp(T2 / t) * t) );
//	printf("function: %lf %lf %lf %lf %lf\n", r, T1, S1, t, ans);
	return ans;
}

double Well_Func_u_deriv ( double t, void *params )
{
	struct u_params *p = (struct u_params *) params;
	double ans;

	double T1 = p->T1;
	double T2 = p->T2;
	double S1 = p->S1;
	double S2 = p->S2;
	double r = p->r;

	ans = ( (pow(r,2) * pow(10,S1) ) / (4 * pow(10,T1)) * exp( ( S2 - T2 ) / t ) * ( ( t + S2 - T2 ) / pow(t,3) ) );
//	printf("deriv: %lf %lf %lf %lf %lf\n", r, T1, S1, t, ans);
	return ans;
	
}

void Well_Func_u_fdf ( double t, void *params, double *y, double *dy)
{

	struct u_params *p = (struct u_params *) params;

	double T1 = p->T1;
	double T2 = p->T2;
	double S1 = p->S1;
	double S2 = p->S2;
	double r = p->r;

	*y = ( 0.03 - ( pow(r,2) * pow(10,S1) * exp(S2 / t) )/( 4 * pow(10,T1) * exp(T2 / t) * t) );
	*dy = ( (pow(r,2) * pow(10,S1) / ( 4 * pow(10,T1) )) * exp( ( S2 - T2 ) / t ) * ( ( t + S2 - T2 ) / pow(t,3) ) );

}

