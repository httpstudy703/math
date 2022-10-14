EXPORT CONST double xcos(double d) {
  double u, s, t = d;
  int ql;

  if (fabsk(d) < TRIGRANGEMAX2) {
    ql = mla(2, rintk(d * M_1_PI - 0.5), 1);
    d = mla(ql, -PI_A2*0.5, d);
    d = mla(ql, -PI_B2*0.5, d);
  } else if (fabsk(d) < TRIGRANGEMAX) {
    double dqh = trunck(d * (M_1_PI / (INT64_C(1) << 23)) - 0.5 * (M_1_PI / (INT64_C(1) << 23)));
    ql = 2*rintk(d * M_1_PI - 0.5 - dqh * (double)(INT64_C(1) << 23))+1;
    dqh *= 1 << 24;

    d = mla(dqh, -PI_A*0.5, d);
    d = mla( ql, -PI_A*0.5, d);
    d = mla(dqh, -PI_B*0.5, d);
    d = mla( ql, -PI_B*0.5, d);
    d = mla(dqh, -PI_C*0.5, d);
    d = mla( ql, -PI_C*0.5, d);
    d = mla(dqh + ql , -PI_D*0.5, d);
  } else {
    ddi_t ddi = rempi(t);
    ql = ((ddi.i & 3) * 2 + (ddi.dd.x > 0) + 7) >> 1;
    if ((ddi.i & 1) == 0) {
      ddi.dd = ddadd2_d2_d2_d2(ddi.dd, dd(mulsign(3.141592653589793116*-0.5, ddi.dd.x > 0 ? 1 : -1),
					  mulsign(1.2246467991473532072e-16*-0.5, ddi.dd.x > 0 ? 1 : -1)));
    }
    d = ddi.dd.x + ddi.dd.y;
    if (xisinf(t) || xisnan(t)) d = SLEEF_NAN;
  }
  
  s = d * d;

  if ((ql & 2) == 0) d = -d;

  double s2 = s * s, s4 = s2 * s2;
  u = POLY8(s, s2, s4,
	    -7.97255955009037868891952e-18,
	    2.81009972710863200091251e-15,
	    -7.64712219118158833288484e-13,
	    1.60590430605664501629054e-10,
	    -2.50521083763502045810755e-08,
	    2.75573192239198747630416e-06,
	    -0.000198412698412696162806809,
	    0.00833333333333332974823815);
  u = mla(u, s, -0.166666666666666657414808);

  u = mla(s, u * d, d);

  return u;
}

EXPORT CONST double xcos_u1(double d) {
  double u;
  Sleef_double2 s, t, x;
  int ql;
  
  d = fabsk(d);

  if (d < TRIGRANGEMAX2) {
    ql = mla(2, rintk(d * M_1_PI - 0.5), 1);
    s = ddadd2_d2_d_d(d, ql * (-PI_A2*0.5));
    s = ddadd_d2_d2_d(s, ql * (-PI_B2*0.5));
  } else if (d < TRIGRANGEMAX) {
    double dqh = trunck(d * (M_1_PI / (INT64_C(1) << 23)) - 0.5 * (M_1_PI / (INT64_C(1) << 23)));
    ql = 2*rintk(d * M_1_PI - 0.5 - dqh * (double)(INT64_C(1) << 23))+1;
    dqh *= 1 << 24;

    u = mla(dqh, -PI_A*0.5, d);
    s = ddadd2_d2_d_d (u,  ql * (-PI_A*0.5));
    s = ddadd2_d2_d2_d(s, dqh * (-PI_B*0.5));
    s = ddadd2_d2_d2_d(s,  ql * (-PI_B*0.5));
    s = ddadd2_d2_d2_d(s, dqh * (-PI_C*0.5));
    s = ddadd2_d2_d2_d(s,  ql * (-PI_C*0.5));
    s = ddadd_d2_d2_d(s, (dqh + ql) * (-PI_D*0.5));
  } else {
    ddi_t ddi = rempi(d);
    ql = ((ddi.i & 3) * 2 + (ddi.dd.x > 0) + 7) >> 1;
    if ((ddi.i & 1) == 0) {
      ddi.dd = ddadd2_d2_d2_d2(ddi.dd, dd(mulsign(3.141592653589793116*-0.5, ddi.dd.x > 0 ? 1 : -1),
					  mulsign(1.2246467991473532072e-16*-0.5, ddi.dd.x > 0 ? 1 : -1)));
    }
    s = ddnormalize_d2_d2(ddi.dd);
    if (xisinf(d) || xisnan(d)) s.x = SLEEF_NAN;
  }
  
  t = s;
  s = ddsqu_d2_d2(s);

  double s2 = s.x * s.x, s4 = s2 * s2;
  u = POLY6(s.x, s2, s4,
	    2.72052416138529567917983e-15,
	    -7.6429259411395447190023e-13,
	    1.60589370117277896211623e-10,
	    -2.5052106814843123359368e-08,
	    2.75573192104428224777379e-06,
	    -0.000198412698412046454654947);
  u = mla(u, s.x, 0.00833333333333318056201922);   //17375093224066086

  x = ddadd_d2_d_d2(1, ddmul_d2_d2_d2(ddadd_d2_d_d(-0.166666666666666657414808, u * s.x), s));
  u = ddmul_d_d2_d2(t, x);
  
  if ((((int)ql) & 2) == 0) u = -u;

  return u;
}
