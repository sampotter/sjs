#include "mat.h"

#include <assert.h>

dmat22 dmat22_add(dmat22 A, dmat22 B) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x + B.rows[0].x, A.rows[0].y + B.rows[0].y},
      {A.rows[1].x + B.rows[1].x, A.rows[1].y + B.rows[1].y},
    }
  };
}

dmat22 dmat22_sub(dmat22 A, dmat22 B) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x - B.rows[0].x, A.rows[0].y - B.rows[0].y},
      {A.rows[1].x - B.rows[1].x, A.rows[1].y - B.rows[1].y},
    }
  };
}

dmat22 dmat22_mul(dmat22 A, dmat22 B) {
  dmat22 C;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      C.data[i][j] = 0;
      for (int k = 0; k < 2; ++k) {
        C.data[i][j] += A.data[i][k]*B.data[k][j];
      }
    }
  }
  return C;
}

dmat22 dmat22_dbl_mul(dmat22 A, dbl a) {
  return (dmat22) {
    .rows = {
      {a*A.rows[0].x, a*A.rows[0].y},
      {a*A.rows[1].x, a*A.rows[1].y}
    }
  };
}

dmat22 dmat22_dbl_div(dmat22 A, dbl a) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x/a, A.rows[0].y/a},
      {A.rows[1].x/a, A.rows[1].y/a}
    }
  };
}

dvec2 dmat22_dvec2_mul(dmat22 A, dvec2 x) {
  return (dvec2) {
    A.rows[0].x*x.x + A.rows[0].y*x.y,
    A.rows[1].x*x.x + A.rows[1].y*x.y
  };
}

dvec2 dmat22_dvec2_solve(dmat22 A, dvec2 b) {
  dbl det = A.data[0][0]*A.data[1][1] - A.data[0][1]*A.data[1][0];
  return (dvec2) {
    .x = (A.data[1][1]*b.x - A.data[0][1]*b.y)/det,
    .y = (A.data[0][0]*b.y - A.data[1][0]*b.x)/det
  };
}

dmat22 dvec2_outer(dvec2 u, dvec2 v) {
  return (dmat22) {
    .rows = {
      {u.x*v.x, u.x*v.y},
      {u.y*v.x, u.y*v.y}
    }
  };
}

void dmat22_invert(dmat22 *A) {
  dbl det = A->data[0][0]*A->data[1][1] - A->data[1][0]*A->data[0][1];
  *A = (dmat22) {
    .rows = {
      { A->rows[1].y/det, -A->rows[0].y/det},
      {-A->rows[1].x/det,  A->rows[0].x/det}
    }
  };
}

dbl dmat22_trace(dmat22 const *A) {
  return A->data[0][0] + A->data[1][1];
}

dbl dmat22_det(dmat22 const *A) {
  return A->data[0][0]*A->data[1][1] - A->data[0][1]*A->data[1][0];
}

dvec2 dmat22_eigvals(dmat22 const *A) {
  dbl tr = dmat22_trace(A);
  dbl disc = tr*tr - 4*dmat22_det(A);
  assert(disc >= 0);
  dbl tmp = sqrt(disc);
  return (dvec2) {(tr + tmp)/2, (tr - tmp)/2};
}

void dmat22_transpose(dmat22 *A) {
  dbl tmp = A->data[1][0];
  A->data[1][0] = A->data[0][1];
  A->data[0][1] = tmp;
}

#define DET33(x00, x01, x02, x10, x11, x12, x20, x21, x22) (    \
  x00*(x11*x22 - x21*x12) -                                     \
  x01*(x10*x22 - x20*x12) +                                     \
  x02*(x10*x21 - x20*x11))

dvec3 dmat33_dvec3_solve(dmat33 A, dvec3 b) {
  dbl a00 = A.data[0][0], a01 = A.data[0][1], a02 = A.data[0][2];
  dbl a10 = A.data[1][0], a11 = A.data[1][1], a12 = A.data[1][2];
  dbl a20 = A.data[2][0], a21 = A.data[2][1], a22 = A.data[2][2];
  dbl det = DET33(a00, a01, a02, a10, a11, a12, a20, a21, a22);
  dbl detx = DET33(b.x, a01, a02, b.y, a11, a12, b.z, a21, a22);
  dbl dety = DET33(a00, b.x, a02, a10, b.y, a12, a20, b.z, a22);
  dbl detz = DET33(a00, a01, b.x, a10, a11, b.y, a20, a21, b.z);
  return (dvec3) {detx/det, dety/det, detz/det};
}

#undef DET33

dvec4 dmat44_dvec4_mul(dmat44 const A, dvec4 const x) {
  dvec4 y;
  for (int i = 0; i < 4; ++i) {
    y.data[i] = dvec4_dot(A.rows[i], x);
  }
  return y;
}

dvec4 dvec4_dmat44_mul(dvec4 const x, dmat44 const A) {
  dvec4 y;
  for (int j = 0; j < 4; ++j) {
    y.data[j] = dvec4_dot(x, dmat44_col(A, j));
  }
  return y;
}

dmat44 dmat44_dmat44_mul(dmat44 const A, dmat44 const B) {
  dmat44 C;
  for (int i = 0; i < 4; ++i) {
    C.rows[i] = dvec4_dmat44_mul(A.rows[i], B);
  }
  return C;
}

dvec4 dmat44_col(dmat44 const A, int j) {
  dvec4 a;
  for (int i = 0; i < 4; ++i) {
    a.data[i] = A.rows[i].data[j];
  }
  return a;
}
