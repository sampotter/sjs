#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

// '1': s(x, y) = 1

dbl s_1(dbl x, dbl y, void *context);
dvec2 grad_s_1(dbl x, dbl y, void *context);

dbl u_1(dbl x, dbl y);
dbl ux_1(dbl x, dbl y);
dbl uy_1(dbl x, dbl y);
dbl uxy_1(dbl x, dbl y);

// 'p': s(x, y) = 1/(1 + 0.133*x - 0.0933*y)

dbl s_p(dbl x, dbl y, void *context);
dvec2 grad_s_p(dbl x, dbl y, void *context);

dbl u_p(dbl x, dbl y);
dbl ux_p(dbl x, dbl y);
dbl uy_p(dbl x, dbl y);
dbl uxy_p(dbl x, dbl y);

// 'v': s(x, y) = 2/(1 + x)

dbl s_v(dbl x, dbl y, void *context);
dvec2 grad_s_v(dbl x, dbl y, void *context);

dbl u_v(dbl x, dbl y);
dbl ux_v(dbl x, dbl y);
dbl uy_v(dbl x, dbl y);
dbl uxy_v(dbl x, dbl y);

// // 'm': s(x, y) = sqrt(sin(x + y)^2 + (x + sin(x + y))^2)

// dbl s_m(dbl x, dbl y, void *context);
// dvec2 grad_s_m(dbl x, dbl y, void *context);

// dbl u_m(dbl x, dbl y);
// dbl ux_m(dbl x, dbl y);
// dbl uy_m(dbl x, dbl y);
// dbl uxy_m(dbl x, dbl y);

// // 'g': s(x, y) = sqrt(4 - 6*y)

// dbl s_g(dbl x, dbl y, void *context);
// dvec2 grad_s_g(dbl x, dbl y, void *context);

// dbl u_g(dbl x, dbl y);
// dbl ux_g(dbl x, dbl y);
// dbl uy_g(dbl x, dbl y);
// dbl uxy_g(dbl x, dbl y);
