/*-----
* 定数・構造体の定義
*-----*/


/*定数*/
/*論理型*/
typedef int bool;
#define TRUE 1
#define FALSE 0

/*WGS-84定数*/
#define PI 3.1415926535898
#define C 299792458e8
#define MUe 3.986005e14 /*地球重力定数[m^3/s^2]*/
#define dOMEGAe 7.2921151467e-5 /*地球自転角速度[rad/s]*/
#define Re 6378137.0 /*地球の長半径[m]*/
#define f (1.0/298.257223563) /*地球の扁平率*/

/*角度の変換*/
#define rad_to_deg(rad) ((rad)*180.0)
#define deg_to_rad(deg) ((deg)/180.0*PI)
#define deg_to_sc(rad) ((rad)/PI)
#define sc_to_deg(sc) ((sc)*PI)

/*取り扱える行列の大きさ*/
#define MAX_N 16 /*観測衛星の上限*/
#define MAX_M 4 /*未知数の最大数*/
#define MAX_PRN 32 /*衛星番号の上限*/

/*時間*/
#define SECONDS_DAY (3600L*24L)
#define SECONDS_WEEK (3600L*24L*7L) 

/*構造体*/

/*時刻を表す構造体*/
typedef struct {
	int week; /*週番号*/
	double sec; /*週初めからの経過時間*/
} wrime;

/*直交座標を表す構造体*/
typedef struct {
	double x;
	double y;
	double z;
} posxyz;
# define SQ(x) ((x)*(x))
# define DIST(a,b) sqrt(SQ(a.x-b.x)+SQ(a.y-b.y)+SQ(a.z-b.z))

/*経緯度を表す構造体*/
typedef struct {
	double lat;
	double lon;
	double hgt;
} posblh;

/*ENU座標を表す構造体*/
typedef struct {
	double e; /*東方向*/
	double n; /*北方向*/
	double u; /*上方向*/
} posenu;
