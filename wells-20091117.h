#define MAXNAME 1024

enum AQUIFER_TYPE {CONFINED = 1, UNCONFINED, LEAKY, LEAKY_UNCONFINED};
enum BOUNDARY_TYPE {NO_BOUNDARY = 1, TYPE_1, TYPE_2};

struct Problem_Info
{
	char file[MAXFILENAME];
	char name[MAXNAME];
};

struct Aquifer_Data
{
	int   co_aquifer;
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
	double *k;
	double *S;
	double *B;
	double *sw;
};

struct Point_Data
{
	int   nP;
	char  **id;
	double *x;
	double *y;
	double *h; 
	double *slope;
	int   *nT;
	double **t;
	double **s;
};

struct Grid_Data
{
	int   nx;		int   ny;
	double xmin;		double xmax;
	double ymin;		double ymax;
};


