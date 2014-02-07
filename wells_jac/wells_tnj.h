#define MAXWELL   50
#define MAXRATE   3000
#define MAXTIME   5000
#define MAXPOINT  100

enum AQUIFER_TYPE {CONFINED = 1, UNCONFINED, LEAKY, LEAKY_UNCONFINED};
enum BOUNDARY_TYPE {NO_BOUNDARY = 1, TYPE_1, TYPE_2};

struct Problem_Info
{
	char file[MAXFILENAME];
	char name[61];
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
	double kobs[MAXPOINT][MAXWELL];
	double kprod[MAXWELL][MAXWELL];
	double Sobs[MAXPOINT][MAXWELL];
	double Sprod[MAXWELL][MAXWELL];
	double B;
};

struct Well_Data
{
	int   nW;
	char  id[MAXWELL][61];
	double x[MAXWELL];
	double y[MAXWELL];
	double r[MAXWELL];
	double h[MAXPOINT];
	int   nQ[MAXWELL];
	double Q[MAXWELL][MAXRATE];
	double t[MAXWELL][MAXRATE];
	double s[MAXWELL][MAXRATE];
	double s_t[MAXRATE];
	double s_j[MAXRATE];
	double m[MAXWELL];
	double k[MAXWELL];
	double S[MAXWELL];
	double B[MAXWELL];
};

struct Point_Data
{
	int   nP;
	char  id[MAXPOINT][61];
	double x[MAXPOINT];
	double y[MAXPOINT];
	double h[MAXPOINT]; 
	double k[MAXPOINT];
	double S[MAXPOINT];
	int   nT[MAXPOINT];
	double t[MAXPOINT][MAXTIME];
	double s[MAXPOINT][MAXTIME];
	double s_t[MAXTIME];
	double s_j[MAXTIME];
	double s_t_abs[MAXTIME];
	double s_j_abs[MAXTIME];
};

struct Grid_Data
{
	int   nx;		int   ny;
	double xmin;		double xmax;
	double ymin;		double ymax;
};

struct Jacob_Info
{
	double t_Jac_obs[MAXPOINT][MAXWELL];
	double t_Jac_prod[MAXWELL][MAXWELL];
	double t_tot_obs[MAXPOINT];
	double t_tot_prod[MAXWELL];
	int n_Jac_obs[MAXPOINT][MAXWELL];
	int n_Jac_prod[MAXWELL][MAXWELL];
};


