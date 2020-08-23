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
  (void) context;
  dbl r = NORM(x, y);
  return (dvec2) {.x = x/r, .y = y/r};
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

dbl uxy_1(dbl x, dbl y) {
  return -ux_1(x, y)*uy_1(x, y)/NORM(x, y);
}

#define VNORM sqrt(VX*VX + VY*VY);
#define VNORMSQ (VX*VX + VY*VY)

#define V0 1.0
#define VX 0.133
#define VY -0.0933

dbl s_p(dbl x, dbl y, void *context) {
  (void) context;
  return 1.0/(V0 + VX*x + VY*y);
}

dbl sx_p(dbl x, dbl y) {
  return -VX*pow(s_p(x, y, NULL), 2.0);
}

dbl sy_p(dbl x, dbl y) {
  return -VY*pow(s_p(x, y, NULL), 2.0);
}

dbl sxy_p(dbl x, dbl y) {
  return 2*VX*VY*pow(s_p(x, y, NULL), 3.0);
}

dvec2 grad_s_p(dbl x, dbl y, void *context) {
  (void) context;
  return (dvec2) {.x = sx_p(x, y), .y = sy_p(x, y)};
}

dbl f(dbl x, dbl y) {
  return 1 + s_p(x, y, NULL)*VNORMSQ*NORMSQ(x, y)/2;
}

dbl fx(dbl x, dbl y) {
  return VNORMSQ*(sx_p(x, y)*NORMSQ(x, y) + 2*s_p(x, y, NULL)*x)/2;
}

dbl fy(dbl x, dbl y) {
  return VNORMSQ*(sy_p(x, y)*NORMSQ(x, y) + 2*s_p(x, y, NULL)*y)/2;
}

dbl fxy(dbl x, dbl y) {
  return VNORMSQ*(sxy_p(x, y)*NORMSQ(x, y) + 2*(sx_p(x, y)*y + sy_p(x, y)*x))/2;
}

dbl dacosh(dbl z) {
  return pow(z - 1, -0.5)*pow(z + 1, -0.5);
}

dbl d2acosh(dbl z) {
  return -z*pow(z - 1, -1.5)*pow(z + 1, -1.5);
}

dbl u_p(dbl x, dbl y) {
  return acosh(f(x, y))/VNORM;
}

dbl ux_p(dbl x, dbl y) {
  return dacosh(f(x, y))*fx(x, y)/VNORM;
}

dbl uy_p(dbl x, dbl y) {
  return dacosh(f(x, y))*fy(x, y)/VNORM;
}

dbl uxy_p(dbl x, dbl y) {
  dbl tmp = f(x, y);
  return (d2acosh(tmp)*fx(x, y)*fy(x, y) + dacosh(tmp)*fxy(x, y))/VNORM;
}

#undef V0
#undef VX
#undef VY

#define V0 0.5
#define VX 0.5
#define VY 0.0

dbl s_v(dbl x, dbl y, void *context) {
  (void) context;
  return 1.0/(V0 + VX*x + VY*y);
}

dbl sx_v(dbl x, dbl y) {
  return -VX*pow(s_v(x, y, NULL), 2.0);
}

dbl sy_v(dbl x, dbl y) {
  return -VY*pow(s_v(x, y, NULL), 2.0);
}

dvec2 grad_s_v(dbl x, dbl y, void *context) {
  (void) context;
  return (dvec2) {.x = sx_v(x, y), .y = sy_v(x, y)};
}

dbl u_v(dbl x, dbl y) {
  return acosh(f(x, y))/VNORM;
}

dbl ux_v(dbl x, dbl y) {
  return dacosh(f(x, y))*fx(x, y)/VNORM;
}

dbl uy_v(dbl x, dbl y) {
  return dacosh(f(x, y))*fy(x, y)/VNORM;
}

dbl uxy_v(dbl x, dbl y) {
  dbl tmp = f(x, y);
  return (d2acosh(tmp)*fx(x, y)*fy(x, y) + dacosh(tmp)*fxy(x, y))/VNORM;
}

#undef V0
#undef VX
#undef VY
