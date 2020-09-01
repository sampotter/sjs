#include "slow.h"

#include <math.h>
#include <stddef.h>

#include "vec.h"

#define NORMSQ(x, y) (x*x + y*y)
#define NORM(x, y) sqrt(NORMSQ(x, y))

dbl s_1(dbl x, dbl y, void *context) {
  (void) x;
  (void) y;
  (void) context;
  return 1.0;
}

dvec2 grad_s_1(dbl x, dbl y, void *context) {
  (void) x;
  (void) y;
  (void) context;
  return dvec2_zero();
}

dbl u_1(dbl x, dbl y) {
  return NORM(x, y);
}

dbl ux_1(dbl x, dbl y) {
  return x/NORM(x, y);
}

dbl uy_1(dbl x, dbl y) {
  return y/NORM(x, y);
}

dbl uxx_1(dbl x, dbl y) {
  dbl ux = ux_1(x, y);
  return (1 - ux*ux)/u_1(x, y);
}

dbl uxy_1(dbl x, dbl y) {
  return -ux_1(x, y)*uy_1(x, y)/NORM(x, y);
}

dbl uyy_1(dbl x, dbl y) {
  dbl uy = uy_1(x, y);
  return (1 - uy*uy)/u_1(x, y);
}

#define VNORM sqrt(VX*VX + VY*VY);
#define VNORMSQ (VX*VX + VY*VY)

#define V0 1.0
#define VX 0.133
#define VY -0.0933

dbl s_p(dbl x, dbl y, void *context) {
  (void)context;
  return 1.0/(V0 + VX*x + VY*y);
}

dbl sx_p(dbl x, dbl y, void *context) {
  (void)context;
  return -VX*pow(s_p(x, y, NULL), 2.0);
}

dbl sy_p(dbl x, dbl y, void *context) {
  (void)context;
  return -VY*pow(s_p(x, y, NULL), 2.0);
}

dbl sxx_p(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VX*VX*pow(s_p(x, y, NULL), 3.0);
}

dbl sxy_p(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VX*VY*pow(s_p(x, y, NULL), 3.0);
}

dbl syy_p(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VY*VY*pow(s_p(x, y, NULL), 3.0);
}

dvec2 grad_s_p(dbl x, dbl y, void *context) {
  return (dvec2) {.x = sx_p(x, y, context), .y = sy_p(x, y, context)};
}

#define SCALE (VNORMSQ/(2*V0))

dbl f_p(dbl x, dbl y) {
  return 1 + SCALE*s_p(x, y, NULL)*NORMSQ(x, y);
}

dbl fx_p(dbl x, dbl y) {
  return SCALE*(sx_p(x, y, NULL)*NORMSQ(x, y) + 2*s_p(x, y, NULL)*x);
}

dbl fy_p(dbl x, dbl y) {
  return SCALE*(sy_p(x, y, NULL)*NORMSQ(x, y) + 2*s_p(x, y, NULL)*y);
}

dbl fxx_p(dbl x, dbl y) {
  return SCALE*(
    sxx_p(x, y, NULL)*NORMSQ(x, y) + 4*sx_p(x, y, NULL)*x + 2*s_p(x, y, NULL));
}

dbl fxy_p(dbl x, dbl y) {
  return SCALE*(
    sxy_p(x, y, NULL)*NORMSQ(x, y) + 2*(sx_p(x, y, NULL)*y + sy_p(x, y, NULL)*x));
}

dbl fyy_p(dbl x, dbl y) {
  return SCALE*(
    syy_p(x, y, NULL)*NORMSQ(x, y) + 4*sy_p(x, y, NULL)*y + 2*s_p(x, y, NULL));
}

#undef SCALE

dbl dacosh(dbl z) {
  return pow(z - 1, -0.5)*pow(z + 1, -0.5);
}

dbl d2acosh(dbl z) {
  return -z*pow(z - 1, -1.5)*pow(z + 1, -1.5);
}

dbl u_p(dbl x, dbl y) {
  return acosh(f_p(x, y))/VNORM;
}

dbl ux_p(dbl x, dbl y) {
  return dacosh(f_p(x, y))*fx_p(x, y)/VNORM;
}

dbl uy_p(dbl x, dbl y) {
  return dacosh(f_p(x, y))*fy_p(x, y)/VNORM;
}

dbl uxx_p(dbl x, dbl y) {
  dbl f = f_p(x, y);
  dbl fx = fx_p(x, y);
  return (d2acosh(f)*fx*fx + dacosh(f)*fxx_p(x, y))/VNORM;
}

dbl uxy_p(dbl x, dbl y) {
  dbl f = f_p(x, y);
  return (d2acosh(f)*fx_p(x, y)*fy_p(x, y) + dacosh(f)*fxy_p(x, y))/VNORM;
}

dbl uyy_p(dbl x, dbl y) {
  dbl f = f_p(x, y);
  dbl fy = fy_p(x, y);
  return (d2acosh(f)*fy*fy + dacosh(f)*fyy_p(x, y))/VNORM;
}

#undef V0
#undef VX
#undef VY

#define V0 0.5
#define VX 0.5
#define VY 0.0

dbl s_v(dbl x, dbl y, void *context) {
  (void)context;
  return 1.0/(V0 + VX*x + VY*y);
}

dbl sx_v(dbl x, dbl y, void *context) {
  (void)context;
  return -VX*pow(s_v(x, y, NULL), 2.0);
}

dbl sy_v(dbl x, dbl y, void *context) {
  (void)context;
  return -VY*pow(s_v(x, y, NULL), 2.0);
}

dbl sxx_v(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VX*VX*pow(s_v(x, y, NULL), 3.0);
}

dbl sxy_v(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VX*VY*pow(s_v(x, y, NULL), 3.0);
}

dbl syy_v(dbl x, dbl y, void *context) {
  (void)context;
  return 2*VY*VY*pow(s_v(x, y, NULL), 3.0);
}

dvec2 grad_s_v(dbl x, dbl y, void *context) {
  return (dvec2) {.x = sx_v(x, y, context), .y = sy_v(x, y, context)};
}

#define SCALE (VNORMSQ/(2*V0))

dbl f_v(dbl x, dbl y) {
  return 1 + SCALE*s_v(x, y, NULL)*NORMSQ(x, y);
}

dbl fx_v(dbl x, dbl y) {
  return SCALE*(sx_v(x, y, NULL)*NORMSQ(x, y) + 2*s_v(x, y, NULL)*x);
}

dbl fy_v(dbl x, dbl y) {
  return SCALE*(sy_v(x, y, NULL)*NORMSQ(x, y) + 2*s_v(x, y, NULL)*y);
}

dbl fxx_v(dbl x, dbl y) {
  return SCALE*(
    sxx_v(x, y, NULL)*NORMSQ(x, y) + 4*sx_v(x, y, NULL)*x + 2*s_v(x, y, NULL));
}

dbl fxy_v(dbl x, dbl y) {
  return SCALE*(
    sxy_v(x, y, NULL)*NORMSQ(x, y) + 2*(sx_v(x, y, NULL)*y + sy_v(x, y, NULL)*x));
}

dbl fyy_v(dbl x, dbl y) {
  return SCALE*(
    syy_v(x, y, NULL)*NORMSQ(x, y) + 4*sy_v(x, y, NULL)*y + 2*s_v(x, y, NULL));
}

#undef SCALE

dbl u_v(dbl x, dbl y) {
  return acosh(f_v(x, y))/VNORM;
}

dbl ux_v(dbl x, dbl y) {
  return dacosh(f_v(x, y))*fx_v(x, y)/VNORM;
}

dbl uy_v(dbl x, dbl y) {
  return dacosh(f_v(x, y))*fy_v(x, y)/VNORM;
}

dbl uxx_v(dbl x, dbl y) {
  dbl f = f_v(x, y);
  dbl fx = fx_v(x, y);
  return (d2acosh(f)*fx*fx + dacosh(f)*fxx_v(x, y))/VNORM;
}

dbl uxy_v(dbl x, dbl y) {
  dbl f = f_v(x, y);
  return (d2acosh(f)*fx_v(x, y)*fy_v(x, y) + dacosh(f)*fxy_v(x, y))/VNORM;
}

dbl uyy_v(dbl x, dbl y) {
  dbl f = f_v(x, y);
  dbl fy = fy_v(x, y);
  return (d2acosh(f)*fy*fy + dacosh(f)*fyy_v(x, y))/VNORM;
}

#undef V0
#undef VX
#undef VY

#undef VNORM
#undef VNORMSQ

dbl s_m(dbl x, dbl y, void *context) {
  (void)context;
  dbl tmp1 = sin(x + y);
  dbl tmp2 = x + tmp1;
  return sqrt(tmp1*tmp1 + tmp2*tmp2);
}

dvec2 grad_s_m(dbl x, dbl y, void *context) {
  dbl tmp1 = sin(2.0*(x + y)) + x*cos(x + y);
  dbl tmp2 = tmp1 + x + sin(x + y);
  dbl s = s_m(x, y, context);
  return (dvec2) {.x = tmp2/s, .y = tmp1/s};
}

dbl u_m(dbl x, dbl y) {
  dbl tmp1 = sin((x + y)/2);
  return x*x/2 + 2*tmp1*tmp1;
}

dbl ux_m(dbl x, dbl y) {
  return x + sin(x + y);
}

dbl uy_m(dbl x, dbl y) {
  return sin(x + y);
}

dbl uxx_m(dbl x, dbl y) {
  return 1 + cos(x + y);
}

dbl uxy_m(dbl x, dbl y) {
  return cos(x + y);
}

dbl uyy_m(dbl x, dbl y) {
  return cos(x + y);
}

#define G0 2.0
#define GX 0.0
#define GY -3.0
#define GNORMSQ (GX*GX + GY*GY)

dbl s_g(dbl x, dbl y, void *context) {
  (void)context;
  return sqrt(G0*G0 + 2*(GX*x + GY*y));
}

dvec2 grad_s_g(dbl x, dbl y, void *context) {
  dbl s = s_g(x, y, context);
  return (dvec2) {.x = GX/s, .y = GY/s};
}

// \bar{S}^2 in Fomel's paper
dbl q_g(dbl x, dbl y) {
  return G0*G0 + GX*x + GY*y;
}

dbl T_g(dbl x, dbl y) {
  dbl q = q_g(x, y);
  return sqrt(q*q - GNORMSQ*NORMSQ(x, y));
}

dbl Tx_g(dbl x, dbl y) {
  dbl q = q_g(x, y);
  dbl T = T_g(x, y);
  return (GX*q - GNORMSQ*x)/T;
}

dbl Ty_g(dbl x, dbl y) {
  dbl q = q_g(x, y);
  dbl T = T_g(x, y);
  return (GY*q - GNORMSQ*y)/T;
}

dbl Txx_g(dbl x, dbl y) {
  dbl Tx = Tx_g(x, y);
  return -(GY*GY + Tx*Tx)/T_g(x, y);
}

dbl Txy_g(dbl x, dbl y) {
  return (GX*GY - Tx_g(x, y)*Ty_g(x, y))/T_g(x, y);
}

dbl Tyy_g(dbl x, dbl y) {
  dbl Ty = Ty_g(x, y);
  return -(GX*GX + Ty*Ty)/T_g(x, y);
}

dbl sig_g(dbl x, dbl y) {
  return sqrt(2*(q_g(x, y) - T_g(x, y))/GNORMSQ);
}

dbl sigx_g(dbl x, dbl y) {
  return (GX - Tx_g(x, y))/(GNORMSQ*sig_g(x, y));
}

dbl sigy_g(dbl x, dbl y) {
  return (GY - Ty_g(x, y))/(GNORMSQ*sig_g(x, y));
}

dbl sigxx_g(dbl x, dbl y) {
  dbl sigx = sigx_g(x, y);
  return -(Txx_g(x, y)/GNORMSQ + sigx*sigx)/sig_g(x, y);
}

dbl sigxy_g(dbl x, dbl y) {
  return -(Txy_g(x, y)/GNORMSQ + sigx_g(x, y)*sigy_g(x, y))/sig_g(x, y);
}

dbl sigyy_g(dbl x, dbl y) {
  dbl sigy = sigy_g(x, y);
  return -(Tyy_g(x, y)/GNORMSQ + sigy*sigy)/sig_g(x, y);
}

dbl u_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sig_sq = sig*sig;
  return sig*q_g(x, y)- GNORMSQ*sig*sig_sq/6;
}

dbl ux_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sig_sq = sig*sig;
  dbl sigx = sigx_g(x, y);
  return GX*sig + q_g(x, y)*sigx - GNORMSQ*sig_sq*sigx/2;
}

dbl uy_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sig_sq = sig*sig;
  dbl sigy = sigy_g(x, y);
  return GY*sig + q_g(x, y)*sigy - GNORMSQ*sig_sq*sigy/2;
}

dbl uxx_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sigx = sigx_g(x, y);
  dbl sigxx = sigxx_g(x, y);
  return 2*GX*sigx + q_g(x, y)*sigxx - GNORMSQ*sig*(sigx*sigx + sig*sigxx/2);
}

dbl uxy_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sigx = sigx_g(x, y);
  dbl sigy = sigy_g(x, y);
  dbl sigxy = sigxy_g(x, y);
  return GX*sigy + GY*sigx + q_g(x, y)*sigxy - GNORMSQ*sig*(sigx*sigy + sig*sigxy/2);
}

dbl uyy_g(dbl x, dbl y) {
  dbl sig = sig_g(x, y);
  dbl sigy = sigy_g(x, y);
  dbl sigyy = sigyy_g(x, y);
  return 2*GY*sigy + q_g(x, y)*sigyy - GNORMSQ*sig*(sigy*sigy + sig*sigyy/2);
}

#undef G0
#undef GX
#undef GY
#undef GNORMSQ
