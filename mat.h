#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "vec.h"

typedef struct {
  union {
    dbl data[2][2];
    dvec2 rows[2];
  };
} dmat22;

dmat22 dmat22_add(dmat22 A, dmat22 B);
dmat22 dmat22_sub(dmat22 A, dmat22 B);
dmat22 dmat22_mul(dmat22 A, dmat22 B);
dmat22 dmat22_dbl_mul(dmat22 A, dbl a);
dmat22 dmat22_dbl_div(dmat22 A, dbl a);
dvec2 dmat22_dvec2_mul(dmat22 A, dvec2 x);
dvec2 dmat22_dvec2_solve(dmat22 A, dvec2 b);
dmat22 dvec2_outer(dvec2 u, dvec2 v);
void dmat22_invert(dmat22 *A);
dbl dmat22_trace(dmat22 const *A);
dbl dmat22_det(dmat22 const *A);
dvec2 dmat22_eigvals(dmat22 const *A);
void dmat22_transpose(dmat22 *A);

typedef struct {
  union {
    dbl data[3][3];
    dvec3 rows[3];
  };
} dmat33;

dvec3 dmat33_dvec3_solve(dmat33 A, dvec3 b);

typedef struct {
  union {
    dbl data[4][4];
    dvec4 rows[4];
  };
} dmat44;

dvec4 dmat44_dvec4_mul(dmat44 const A, dvec4 const x);
dvec4 dvec4_dmat44_mul(dvec4 const x, dmat44 const A);
dmat44 dmat44_dmat44_mul(dmat44 const A, dmat44 const B);
dvec4 dmat44_col(dmat44 const A, int j);

#ifdef __cplusplus
}
#endif
