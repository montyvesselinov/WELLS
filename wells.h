#define MAXNAME 1024
#define MAXFILENAME 1024 

enum AQUIFER_TYPE {CONFINED = 1, UNCONFINED, LEAKY, LEAKY_UNCONFINED};
enum BOUNDARY_TYPE {NO_BOUNDARY = 1, TYPE_1, TYPE_2};
/*Solution Type */
enum SOLUTION_TYPE {THEIS = 1, HANTUSH = 2, PAPADOPULOS = 3, Mishra_Neuman = 4, THEIS_OLD = 5};

struct Problem_Info
{
	char file[MAXFILENAME];
	char name[MAXNAME];
};

struct Aquifer_Data
{
	int   co_aquifer;
/* Solution Type */
	int   sol_type;
	int   co_x_axis;
	int   co_y_axis;
	double origin_x;
	double origin_y;
	double origin_angle;
	double h;    
	double m;
	double **kobs;
	double **kprod;
	double **Sobs;
	double **Sprod;
	double B;
};

struct Well_Data
{
	int   nW;
	int   *pType;
	char  **id;
	double *x;
	double *y;
	double *r;
	double *h;
	int   *nQ;
	double **Q;
	double **t;
	double **s;
	double *m;
	double *l_unsat;
	double *k;
	double *S;
	double *kD;
	double *B;
	double *sw;
/*Partial Penetration Variable */
        double *d;
        double *l;
/*Wellbore Storage Parameter */
        double *CwD;       
/*Unsaturated Parameters */
        double *sy;
	double *akD;
	double *acD;
	double *psiD;
};

struct Point_Data
{
	int   nP;
	char  **id;
	double *x;
	double *y;
	/*Z-coordinate below the Top*/
	double *z1;
    double *z2; 
	double *h; 
	double *slope;
	int   *nT;
	double **t;
	double **s;
    int *nZeros;
};

struct Grid_Data
{
	int   nx;    int   ny;
	double xmin; double xmax;
	double ymin; double ymax;
};

struct hankel_data
{
  int nZeros_min; 
  int nZeros_max;
  double EPSABS ;
  double EPSREL;
  int NUMITER;
};

struct laplace_data
{
   double DHTOL, DHALPHA, time_max ;
};


struct params 
{ 
double kD;
double lD; double dD;
double rwD;double cwD;double beta2;
double akD; double acD;double ts;
double ts1;double ts2;double Q1;double Q2;
double l_unsat;
double sigma;double psiD;
double zD1; double zD2;
double DHALPHA, DHTOL, time_max;
};
