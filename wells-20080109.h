#define MAXWELL   50
#define MAXRATE   2000
#define MAXTIME   2000
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
	float origin_x;
	float origin_y;
	float origin_angle;
	float h;
	float m;
	float k;
	float S;
	float B;
};

struct Well_Data
{
	int   nW;
	char  id[MAXWELL][61];
	float x[MAXWELL];
	float y[MAXWELL];
	float r[MAXWELL];
	int   nQ[MAXWELL];
	float Q[MAXWELL][MAXRATE];
	float t[MAXWELL][MAXRATE];
	float s[MAXWELL][MAXRATE];
};

struct Point_Data
{
	int   nP;
	char  id[MAXPOINT][61];
	float x[MAXPOINT];
	float y[MAXPOINT];
	int   nT[MAXPOINT];
	float t[MAXPOINT][MAXTIME];
	float s[MAXPOINT][MAXTIME];
};

struct Grid_Data
{
	int   nx;		int   ny;
	float xmin;		float xmax;
	float ymin;		float ymax;
};


