/*---------------------------------------
* 航法メッセージの処理
* --------------------------------------*/
#include "constants.h"
/* RINEXファイルの情報*/
#define RINEX_POS_COMMENT 60
#define RINEX_NAV_LINES 8
#define RINEX_NAV_FIELDS_LINE 4

/* 記憶するエフェメリスの最大数 */
#define MAX_EPHMS 20

/* エフェメリスの有効期限[h] */
#define EPHMERIS_EXPIRE 2.0

/* エフェメリスを格納するための構造体 */
typedef struct {
	int week; /* 週番号 */
	double data[RINEX_NAV_FIELDS_LINE * RINEX_NAV_LINES]; /* エフェメリスデータ */
} ephm_info;

/* パラメータ番号を定義 */
/*
* line1 :
			TOC: 時計データの基準時刻
			AF0, AF1, AF2: 衛星時計補正係数（定数項、1次項、2次項）
* line2:
			IODE: エフェメリスデータの発行番号
      Crs: 軌道半径の正弦補正項
      dn: 平均運動差
      M0: 平均近点角 
* line3:
			Cuc, Cus: 緯度引数の余弦・正弦補正項
      e: 離心率
      sqrtA: 軌道長半径の平方根
* line4:
			TOE: エフェメリスの基準時刻
      Cic, Cis: 軌道傾斜角の余弦・正弦補正項 
* line5:
			i0: 軌道面の傾斜角
			Crc: 軌道半径の余弦補正項
			omega: 昇交点赤経
			dOmega: 昇交点赤経の変化率
* line6:
			di: 軌道面の変化率
			CAonL2: L2のコード位相
			WEEK: 週番号
			L2P: L2の精密コード位相
* line7:
			acc: 衛星クロック補正
			health: 衛星の健康状態
			TGD: グループ遅延
			IODC: クロック補正の発行番号
* line8:
			TOT: 送信時刻
			FIT: データ有効期限
*/
enum ephm_para {
	EPHM_TOC, EPHM_AFO, EPHM_AF1, EPHM_AF2, /* line 1 */
	EPHM_IODE, EPHM_Crs, EPHM_d_n, EPHM_MO, /* line 2 */
	EPHM_Cuc, EPHM_e, EPHM_Cus, EPHM_sqrtA, /* line 3 */
	EPHM_TOE, EPHM_Cic, EPHM_OMEGAO, EPHM_Cis, /* line 4 */
	EPHM_i0, EPHM_CrC, EPHM_omega, EPHM_dOmega, /* line 5 */
	EPHM_di, EPHM_CAonL2, EPHM_WEEK, EPHM_L2P, /* line 6 */
	EPHM_acc, EPHM_health, EPHM_TGD, EPHM_IODC, /* line 7 */
	EPHM_TOT, EPHM_FIT/* line 8 */
	};

/* エフェメリスを保持する配列 */
/*
 * ephm_buf: 各衛星（MAX_PRN）の複数エフェメリス（MAX_EPHMS）を格納
 * ephm_count: 各衛星のエフェメリス数を記録
 * current_ephm: 各衛星の現在使用中のエフェメリスのインデックス
 * current_week: 現在のGPS週を記録（初期値-1は未設定を示す）
*/
static ephm_info ephm_buf[MAX_PRN][MAX_EPHMS];
static int ephm_count[MAX_PRN];
static int current_ephm[MAX_PRN];
static int current_week=-1;

/* 閏秒の情報 */
static int leap_sec = 0;

/* 電離層補正情報 */
#define IONO_PARAMETERS 4
static double iono_alpha[IONO_PARAMETERS];
static double iono_beta[IONO_PARAMETERS];

/* 文字列⇨実装の変換 */
static double atof2(char* str)
{
	char* p;

	/* 'D'を'E'に変換する */
	for (p = str; *p != '\0'; p++) {
		if (*p == 'D' || *p == 'd') *p = 'E';
	}

	/* 実数に変換する */
	return atof(str);
}

/* コメント情報を調べる
 * is_comment() - RINEXファイルの60行め以降にあるラベルを調べる
 */
static bool is_comment(char* str)
{
	return (strncmp(linebuf+RINEX_POS_COMMENT, str, strlen(str)) == 0);
}

/*------------------------------------
* read_RINEX_NAV() RINEX航法ファイルを読み込む
--------------------------------------
* void read_RINEX_NAV(fp):
* 	FILE *fp; RINEX航法ファイルのポインタ
--------------------------------------*/

void read_RINEX_NAV(FILE* fp)
{
	int i,j,n,prn,line;
	bool noerr=FALSE;
	double d;
	wtime wt;
	struct tm tmbuf;
	ephm_info info;

	/* 初期化 */
	if (current_week < 0) {
		for (i=0; i<IONO_PARAMETERS; i++) {
			iono_alpha[i] = 0.0;
			iono_beta[i] = 0.0;
		}
		for (i=0; i<MAX_PRN; i++) {
			ephm_count[i] = 0;
		}
	}

	/* ヘッダ部分 */
	while (read_line(fp)) {
		if (is_comment("ION ALPHA")) {
			iono_alpha[0] = atof(get_field(14));
			iono_alpha[1] = atof(get_field(12));
			iono_alpha[2] = atof(get_field(12));
			iono_alpha[3] = atof(get_field(12));
		} else if (is_comment("ION BETA")) {
			iono_beta[0] = atof(get_field(14));
			iono_beta[1] = atof(get_field(12));
			iono_beta[2] = atof(get_field(12));
			iono_beta[3] = atof(get_field(12));
		} else if (is_comment("LEAP SECONDS")) {
			leap_sec = atoi(get_field(6));
		} else if (is_comment("END OF HEADER")) {
			break;
		}
	}

	/* 本文を読む */
	fprintf(stderr, "Reading RINEX NAV ...");
	while(read_line(fp)) {
		n=0;
		for (line=0; line<RINEX_NAV_LINES; line++) {
			/* 左端のデータを読む */
			if (line==0) {
				/* 最初の行 */
				prn = atoi(get_field(2));
				tmbuf.tm_year = atoi(get_field(3));
				if (tmbuf.tm_year < 80) tmbuf.tm_year += 100;
				tmbuf.tm_mon = atoi(get_field(3))-1;
				tmbuf.tm_mday = atoi(get_field(3));
				tmbuf.tm_hour = atoi(get_field(3));
				tmbuf.tm_min = atoi(get_field(3));
				tmbuf.tm_sec = 0;
				wt = date_to_wtime(tmbuf);
				wt.sec += atof(get_field(5));
				info.week = wt.week;
				d = wt.sec;
			} else {
				/* 2行目以降 */
				if (!read_line(fp)) goto err;
				d = atof(get_field(22));
			}
			info.data[n++] = d;

			/* 残りのデータを読む */
			for (i=1, i<RINEX_NAV_FIELDS_LINE; i++) {
				d = atof(get_field(19));
				info.data[n++] = d;
		}
	}
	if (prn<1) continue;

	/* 週番号を揃える */
	t = (info.week-info.data[EPHM_WEEK])*SECONDS_WEEK;
	info.week = info.data[EPHM_WEEK];
	info.data[EPHM_TOE] += t;
	current_week = info.week;

	/* すでに同じものがないかどうか */
	for (i=0; i<ephm_count[prn-1]; i++) {
		/* 週番号が一致するものが対象 */
		if (ephm_buf[prn-1][i].week == info.week) continue;

		/* IODCが一致しているか */
		if (ephm_buf[prn-1][i].data[EPHM_IODC] == info.data[EPHM_IODC]) {
			/* 送信時刻の早いものを残す */
			if (info.data[EPHM_TOT] < ephm_buf[prn-1][i].data[EPHM_TOT]) {
				ephm_buf[prn-1][i] = info;
			}
			prn=0;
			break;
		}
	}
	if (prn<1) continue;

	/* 配列に格納する */
	if (ephm_count[prn-1] >= MAX_EPHMS) {
		fprintf(stderr, "Too Long Nav File. \n");
		goto NOERROR;
	}
	for (i=0; i<ephm_count[prn-1]; i++) {
		t = (ephm_buf[prn-1][i].week-info.week)*SECONDS_WEEK + ephm_buf[prn-1][i].data[EPHM_TOT]
		if (info.data[EPHM_TOT] < t) break;
	}
	for (j=ephm_count[prn-1]; j>i; j--) {
		ephm_buf[prn-1][j] = ephm_buf[prn-1][j-1];
	}
	ephm_buf[prn-1][i] = info;
	ephm_count[prn-1]++;
}
NOERROR:
	noerr = TRUE;

ERROR:
	/* エフェメリスが無い場合 */
	if (current_week < 0) {
		fprintf(stderr, "No ephemeris information. \n");
		return;
	}

	/* エフェメリスのある衛星の数 */
	n=0; for(prn=1; prn<=MAX_PRN; prn++) {
		if (ephm_count[prn-1] > 0) n++;
	}
	fprintf(stderr, "week %d: %d satellites\n", info.week, n);

	/* 途中でファイルが終わっていた場合 */
	if (!noerr) {
		fprintf(stderr, "Error: Unexpected EOF. \n");
	}
