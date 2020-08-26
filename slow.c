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

dbl S_g(dbl x, dbl y) {
  return sqrt(G0*G0 + GX*x + GY*y);
}

dbl Sx_g(dbl x, dbl y) {
  return GX/(2*S_g(x, y));
}

dbl Sy_g(dbl x, dbl y) {
  return GY/(2*S_g(x, y));
}

dbl Sxx_g(dbl x, dbl y) {
  dbl Sx = Sx_g(x, y);
  return -Sx*Sx/S_g(x, y);
}

dbl Sxy_g(dbl x, dbl y) {
  return -Sx_g(x, y)*Sy_g(x, y)/S_g(x, y);
}

dbl Syy_g(dbl x, dbl y) {
  dbl Sy = Sy_g(x, y);
  return -Sy*Sy/S_g(x, y);
}

dbl T_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl S_sq = S*S;
  return sqrt(S_sq*S_sq - GNORMSQ*(x*x + y*y));
}

dbl Tx_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  return (GX*S*S - GNORMSQ*x)/T_g(x, y);
}

dbl Ty_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  return (GY*S*S - GNORMSQ*y)/T_g(x, y);
}

dbl Txx_g(dbl x, dbl y) {
  dbl T = T_g(x, y);
  dbl Tx = Tx_g(x, y);
  return (GX*GX - GNORMSQ - Tx*Tx)/T;
}

dbl Txy_g(dbl x, dbl y) {
  return (GX*GY - Tx_g(x, y)*Ty_g(x, y))/T_g(x, y);
}

dbl Tyy_g(dbl x, dbl y) {
  dbl T = T_g(x, y);
  dbl Ty = Ty_g(x, y);
  return (GY*GY - GNORMSQ - Ty*Ty)/T;
}

dbl sigma_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl T = T_g(x, y);
  return sqrt(2*(S*S - T))/sqrt(GNORMSQ);
}

dbl sigmax_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl S = S_g(x, y);
  dbl T = T_g(x, y);
  dbl Tx = Tx_g(x, y);
  return (x/sig - sig*(GX + Tx)/2)/(S*S + T);
}

dbl sigmay_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl S = S_g(x, y);
  dbl T = T_g(x, y);
  dbl Ty = Ty_g(x, y);
  return (y/sig - sig*(GY + Ty)/2)/(S*S + T);
}

dbl sigmaxx_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl sigx = sigmax_g(x, y);
  dbl T = T_g(x, y);
  dbl Tx = Tx_g(x, y);
  dbl Txx = Txx_g(x, y);
  dbl S = S_g(x, y);
  return -(sigx*sigx/sig + (2*sigx*(GX + Tx) + sig*Txx/2)/(S*S + T));
}

dbl sigmaxy_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl sigx = sigmax_g(x, y);
  dbl sigy = sigmay_g(x, y);
  dbl T = T_g(x, y);
  dbl Tx = Tx_g(x, y);
  dbl Ty = Ty_g(x, y);
  dbl Txy = Txy_g(x, y);
  dbl S = S_g(x, y);
  return -(sigx*sigy/sig + (sigx*(GY + Ty) + sigy*(GX + Tx) + sig*Txy/2)/(S*S + T));
}

dbl sigmayy_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl sigy = sigmay_g(x, y);
  dbl T = T_g(x, y);
  dbl Ty = Ty_g(x, y);
  dbl Tyy = Tyy_g(x, y);
  dbl S = S_g(x, y);
  return -(sigy*sigy/sig + (2*sigy*(GY + Ty) + sig*Tyy/2)/(S*S + T));
}

dbl u_g(dbl x, dbl y) {
  dbl sig = sigma_g(x, y);
  dbl S = S_g(x, y);
  return sig*(S - GNORMSQ*sig*sig)/6;
}

dbl ux_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl Sx = Sx_g(x, y);
  dbl sig = sigma_g(x, y);
  dbl sigx = sigmax_g(x, y);
  return S*sigx + Sx*sig - GNORMSQ*sigx*sig*sig/2;
}

dbl uy_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl Sy = Sy_g(x, y);
  dbl sig = sigma_g(x, y);
  dbl sigy = sigmay_g(x, y);
  return S*sigy + Sy*sig - GNORMSQ*sigy*sig*sig/2;
}

dbl uxx_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl Sx = Sx_g(x, y);
  dbl Sxx = Sxx_g(x, y);
  dbl sig = sigma_g(x, y);
  dbl sigx = sigmax_g(x, y);
  dbl sigxx = sigmaxx_g(x, y);
  return S*sigxx + 2*Sx*sigx + Sxx*sig - GNORMSQ*(sigxx*sig*sig/2 + sig*sigx);
}

dbl uxy_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl Sx = Sx_g(x, y);
  dbl Sy = Sy_g(x, y);
  dbl Sxy = Sxy_g(x, y);
  dbl sig = sigma_g(x, y);
  dbl sigx = sigmax_g(x, y);
  dbl sigy = sigmay_g(x, y);
  dbl sigxy = sigmaxy_g(x, y);
  return S*sigxy + Sy*sigx + Sx*sigy + Sxy*sig
    - GNORMSQ*(sigxy*sig*sig + 2*sigx*sigy*sig)/2;
}

dbl uyy_g(dbl x, dbl y) {
  dbl S = S_g(x, y);
  dbl Sy = Sy_g(x, y);
  dbl Syy = Syy_g(x, y);
  dbl sig = sigma_g(x, y);
  dbl sigy = sigmay_g(x, y);
  dbl sigyy = sigmayy_g(x, y);
  return S*sigyy + 2*Sy*sigy + Syy*sig - GNORMSQ*(sigyy*sig*sig/2 + sig*sigy);
}

#undef G0
#undef GX
#undef GY
#undef GNORMSQ
