#define MAXNAME 1024
#define MAXFILENAME  1024
#define TRUE    1
#define FALSE   0
#define YES     1
#define NO      0
#define READ        0
#define WRITE       1
#define APPEND      2

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
	double **k2obs;
	double **kprod;
	double **Sobs;
	double **S2obs;
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
	double *k2;
	double *S;
	double *S2;
	double *B;
	double *sw;
	int	*kfcn;
	int	*Sfcn;
};

struct Point_Data
{
	int   nP;
	char  **id;
	double *x;
	double *y;
	double *h; 
	double *c0;
	double *c1;
	int   *nT;
	double **t;
	double **s;
	int   trnd;
};

struct Grid_Data
{
	int   nx;		int   ny;
	double xmin;		double xmax;
	double ymin;		double ymax;
};


