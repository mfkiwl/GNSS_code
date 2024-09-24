#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "constants.h"

/*-----------
* 測位計算を行うコード
-----------*/


/*
* 逆行列の計算をするプログラム
* invert_matrix() - 逆行列を計算する
* inverse_matrix(a,n):
* 	double a[][];元の行列(逆行列で上書きされる)
* 	int m;行列の次数
*/

static void inverse_matrix(double a[MAX_M][MAX_M], int m)
{
	int i,j,k;
	double b[MAX_M][MAX_M + MAX_M];

	/*操作用の行列を作成する*/
	for(i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			b[i][j] = a[i][j];
			b[i][j + m] = (i == j) ? 1.0 : 0.0;
		}
	}

	/*ガウスの消去法*/
	for(i = 0; i < m; i++) {
		/*第一行をb[i][j]で正規化する*/
		if (fabs(b[i][i]) <=1E-10) {
			fprintf(stderr, "inverse_matrix: singular matrix\n");
			exit(2);
		}
		for (j= m+m-1; j >= i; j--) {
			b[i][j] /= b[i][i];
		}

		/*他の行の第i列を消去する*/
		for (k=0; k < m; k++) {
			if (k != i) {
				for (j = m+m-1; j >= i; j--) {
					b[k][j] -= b[k][i] * b[i][j];
				}
			}
		}

		/*元の行列を逆行列で上書きする*/
		for (i=0; i < m; i++) {
			for (j=0; j < m; j++) {
				a[i][j] = b[i][j+m];
			}
		}
	}
}

/*
* computer_solution() - 最小二乗法で解を求める
* void computer_solution(G, dr, wgt, dx, cov, n, m):
* 	double G[][];デザイン行列
* 	double dr[];方程式の右辺(n次元)
* 	double wgt[];重み(n次元/NULLなら重みはない)
* 	double dx[];解(m次元)
* 	double cov[][];共分散行列(m*m次元)
* 	int n;方程式の数
* 	int m;未知数の数
* 	
* 与えられた方程式を最小二乗法で解く。重みが不要な場合、はwgtにNULLを指定する。
*/

void computer_solution(double G[MAX_N][MAX_M], double dr[MAX_N], double wgt[MAX_N], double dx[MAX_M], double cov[MAX_M][MAX_M], int n, int m)
{
	int i,j,k;
	double w,a[MAX_M][MAX_M];

	/*GtG を求める*/
	for (i=0; i < m; i++) {
		for (j=0; j < m; j++) {
			cov[i][j] = 0.0;
			for (k=0; k < n; k++) {
				if (wgt==NULL) w=1.0; else w = wgt[k];
				cov[i][j] += G[k][i] * G[k][j] * w;
			}
		}
	}

	/* 逆行列を求める(これが共分散行列Cになる)*/
	inverse_matrix(cov, m);

	/*Gt をかける*/
	for (i=0; i < m; i++) {
		for (j=0; j < n; j++) {
			a[i][j] = 0.0;
			for (k=0; k < m; k++) {
				if (wgt==NULL) w=1.0; else w = wgt[j];
				a[i][j] += cov[i][k] * G[j][k] * w;
			}
		}
	}

	/*drをかけると解になる*/
	for (i=0; i < m; i++) {
		dx[i] = 0.0;
		for (j=0; j < n; j++) {
			dx[i] += a[i][j] * dr[j];
		}
	}
}

#define LOOP 8
#define SATS 5

static posxyz positon[SATS] = {
	{ -13897607.6294, -10930188.6233, 19676689.6804}, /* PRN 05 */
	{ -17800899.1998, 15689920.8120, 11943543.3888}, /* PRN 14 */
	{ -1510958.2282, 26280096.7818,  -3117646.1949}, /* PRN 16 */
	{ -12210758.3517, 20413597.0201, -11649499.5474}, /* PRN 22 */
	{ -170032.6981, 17261822.6784, 20555984.4061},  /* PRN 25 */
};

/* 受信機からの距離 */
static double range[SATS] = {
	23634878.5219, /* PRN 05 */
	20292688.3557, /* PRN 14 */
	24032055.0372, /* PRN 16 */
	24383229.3740, /* PRN 22 */
	22170992.8187, /* PRN 25 */
};

void main(int argc, char **argv)
{
	int i,n,loop;
	double r, G[MAX_N][MAX_M], dr[MAX_N], dx[MAX_M];
	double sol[MAX_M], cov[MAX_M][MAX_M];
	posxyz satpos;

	/* 解を初期化 */
	for (i=0; i<MAX_M; i++) {
		sol[i] = 0.0;
	}

	/* 解を求めるループ*/
	for(loop=0, loop<LOOP; loop++) {
		n = SATS;
		for(i=0; i<n; i++) {
			satpos = positon[i];
			
			/* デザイン行列を作成する */
			r = sqrt((satpos.x - sol[0]) * (satpos.x - sol[0]) +
					 (satpos.y - sol[1]) * (satpos.y - sol[1]) +
					 (satpos.z - sol[2]) * (satpos.z - sol[2]));
			G[i][0] = (sol[0] - satpos.x) / r;
			G[i][1] = (sol[1] - satpos.y) / r;
			G[i][2] = (sol[2] - satpos.z) / r;

			/* 擬似距離の修正量 */
			dr[i] = range[i] - r;
		}

		/* 最小二乗法で解を求める */
		computer_solution(G, dr, NULL, dx, cov, n, 3);

		/* 初期値を加える */
		for (i=0; i<3; i++) {
			sol[i] += dx[i];
		}
		/* 途中経過を出力する */
		printf("LOOP %d: X = %.4f, Y = %.4f, Z = %.4f\n", loop+1, sol[0], sol[1], sol[2]);
	}
	exit(0);
}
