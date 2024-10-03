/*---------------------------------------
* set_ephemeris() - エフェメリスをセット
* --------------------------------------
* bool set_ephemeris(prn, wt, iode); TRUE:セットした
* int prn; 衛星PRN番号 (1~)
* wtime wt; 時刻を指定
* int iode; IODEを指定/ -1:指定なし
* --------------------------------------*/

bool set_ephemeris(int prn, wtime wt, int iode)
{
	int i;
	double t, t0;

	/* 新しいエフェメリスから探す */
	for (i=ephm_count[prn-1]; i>=0; i--) {
		/* 有効期限内であること */
		t0 = (ephm_buf[prn-1][i].week - current_week) * SECONDS_WEEK;
		t = ephm_buf[prn-1][i].data[EPHM_TOC] + t0;
		if (wt.sec < t - EPHEMERIS_EXPIRE*3600.0 + 0.1 < wt.sec) continue;

		/* IODEをチェック */
		if (iode >= 0) {
			if (ephm_buf[prn-1][i].data[EPHM_IODE] == iode) break;
			continue;
		}

		/* 受信時刻 */
		t = ephm_buf[prn-1][i].data[EPHM_TOC] + t0;
		if (t < wt.sec + 0.1) break; /* wtよりも前の時刻 */
	}

	/* カレント情報としてセットする */
	current_ephm[prn-1] = i;
	return (i >= 0);
}

/*---------------------------------------
* get_ephemeris() - エフェメリスのパラメータを取得
* --------------------------------------
* double get_ephemeris(prn, para); パラメータ値
* int prn; 衛星PRN番号 (1~)
* int para; パラメータ番号 (0~)
* --------------------------------------
* 事前にset_ephemeris()でエフェメリスがセットされておくこと
* --------------------------------------*/

double get_ephemeris(int prn, int para)
{
	if (ephm_count[prn-1] < 1 || current_ephm[prn-1]) {
		fprintf(stderr, "Missing ephemeris: PRN=%d.\n", prn);
		exit(2);
	}

	return ephm_buf[prn-1][current_ephm[prn-1]].data[para];
}

/*---------------------------------------
* satellite_clock() - 衛星クロック誤差を計算(不完全版?)
* --------------------------------------
* double satellite_clock(prn, wt); 衛星クロック誤差[s]
* int prn; 衛星PRN番号 (1~)
* wtime wt; 時刻を指定
* --------------------------------------
* 事前にset_ephemeris()でエフェメリスがセットされておくこと
* 内部でget_ephemeris()を呼び出す
* --------------------------------------*/

double satellite_clock(int prn, wtime wt)
{
	double tk, tk0, dt, tr=0.0;

	/* 時間差を求める */
	tk0 = (wt.week - get_ephemeris(prn, EPHM_WEEK)) * SECONDS_WEEK + wt.sec - get_ephemeris(prn, EPHM_TOC);
	tk = tk0; /* 後で利用 */

	/* 衛星時計の補正量を計算 */
	dt = get_ephemeris(prn, EPHM_AF0) + get_ephemeris(prn, EPHM_AF1) * tk + get_ephemeris(prn, EPHM_AF2) * tk * tk;

	return dt + tr - get_ephemeris(prn, EPHM_TGD);
}

/*---------------------------------------
* satellite_position() - 衛星位置を計算(不完全版?)
* --------------------------------------
* posxyz satellite_position(prn, wt); 衛星位置
* int prn; 衛星PRN番号 (1~)
* wtime wt; 時刻を指定
* --------------------------------------
* 事前にset_ephemeris()でエフェメリスがセットされておくこと
* 内部でget_ephemeris()を呼び出す
* --------------------------------------*/

posxyz satellite_position(int prn, wtime wt)
{
	int i;
	double tk, tk0, sqrtA, e, n, Ek, Mk, xk, yk, Omegak,
		vk, pk, uk, rk, ik, d_uk, d_rk, d_ik;
	posxyz pos;

	/* 時間差を求める */
	tk0 = (wt.week - get_ephemeris(prn, EPHM_WEEK)) * SECONDS_WEEK + wt.sec - get_ephemeris(prn, EPHM_TOE);
	tk = tk0; /* 後で利用 */

	/* 離心近点角 Ek[rad] を求める */
	sqrtA = sqrt(get_ephemeris(prn, EPHM_sqrtA));
	e = get_ephemeris(prn, EPHM_e); /* 離心率 */
	n = sqrt(MUe)/sqrtA/sqrtA/sqrtA+get_ephemeris(prn, EPHM_d_n); /* 角速度 */
	Mk = get_ephemeris(prn, EPHM_M0) + n * tk; /* 平均近点角 Mk[rad] */
	Ek = Mk; for(i=0; i<10; i++) Ek = Mk + e * sin(Ek); /* Kepler方程式 */

	/* 軌道面内における衛星位置 */
	rk = sqrtA * sqrt(1.0 - e * cos(Ek)); /* 動径長 */
	vk = atan2((sqrt(1.0 - e * e) * sin(Ek)), (cos(Ek) - e)); /* 真近点角 */
	pk = vk + get_ephemeris(prn, EPHM_omega); /* 緯度引数[rad] */

	/* 補正係数を適用する */
	d_uk = get_ephemeris(prn, EPHM_Cuc) * cos(2.0 * pk) + get_ephemeris(prn, EPHM_Cus) * sin(2.0 * pk);
	d_rk = get_ephemeris(prn, EPHM_Crc) * cos(2.0 * pk) + get_ephemeris(prn, EPHM_Crs) * sin(2.0 * pk);
	d_ik = get_ephemeris(prn, EPHM_Cic) * cos(2.0 * pk) + get_ephemeris(prn, EPHM_Cis) * sin(2.0 * pk);

	uk = pk + d_uk; /* 緯度引数 [rad] */
	rk = rk + d_rk; /* 動径長 [m] */
	ik = get_ephemeris(prn, EPHM_i0) + d_ik + get_ephemeris(prn, EPHM_di) * tk; /* 軌道面の傾斜角 [rad] */

	/* 軌道面内での位置 */
	xk = rk * cos(uk);
	yk = rk * sin(uk);

	/* 昇交点の緯度 [rad] */
	Omegak = get_ephemeris(prn, EPHM_OMEGAO) + (get_ephemeris(prn, EPHM_dOmega) - dOMEGAe) * tk0 - dOMEGAe * get_ephemeris(prn, EPHM_TOE);

	/* ECEF 座標系に変換 */
	pos.x = xk * cos(Omegak) - yk * cos(ik) * sin(Omegak);
	pos.y = xk * sin(Omegak) + yk * cos(ik) * cos(Omegak);
	pos.z = yk * sin(ik);

	return pos;
}