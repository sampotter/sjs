#include "eik.h"
#include "npy.h"
#include "slow.h"

#include <stdlib.h>
#include <stdio.h>

#define MAX(x, y) x > y ? x : y

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("usage: %s <slo_fun> <N>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char slo_fun = argv[1][0];

  dbl (*s)(dbl, dbl, void *);
  dvec2 (*grad_s)(dbl, dbl, void *);

  dbl (*u)(dbl, dbl);
  dbl (*ux)(dbl, dbl);
  dbl (*uy)(dbl, dbl);
  dbl (*uxx)(dbl, dbl);
  dbl (*uxy)(dbl, dbl);
  dbl (*uyy)(dbl, dbl);

  if (slo_fun == '1') {
    s = s_1;
    grad_s = grad_s_1;
    u = u_1;
    ux = ux_1;
    uy = uy_1;
    uxx = uxx_1;
    uxy = uxy_1;
    uyy = uyy_1;
  } else if (slo_fun == 'p') {
    s = s_p;
    grad_s = grad_s_p;
    u = u_p;
    ux = ux_p;
    uy = uy_p;
    uxx = uxx_p;
    uxy = uxy_p;
    uyy = uyy_p;
  } else if (slo_fun == 'v') {
    s = s_v;
    grad_s = grad_s_v;
    u = u_v;
    ux = ux_v;
    uy = uy_v;
    uxx = uxx_v;
    uxy = uxy_v;
    uyy = uyy_v;
  } else if (slo_fun == 'm') {
    s = s_m;
    grad_s = grad_s_m;
    u = u_m;
    ux = ux_m;
    uy = uy_m;
    uxx = uxx_m;
    uxy = uxy_m;
    uyy = uyy_m;
  } else if (slo_fun == 'g') {
    s = s_g;
    grad_s = grad_s_g;
    u = u_g;
    ux = ux_g;
    uy = uy_g;
    uxx = uxx_g;
    uxy = uxy_g;
    uyy = uyy_g;
  } else {
    fprintf(
      stderr,
      "ERROR: slo_fun should be one of: 1, p, v, m, g (got '%c')\n",
      slo_fun
      );
    exit(EXIT_FAILURE);
  }

  eik_s * scheme;
  eik_alloc(&scheme);

  field2_s slow = {
    .f = s,
    .grad_f = grad_s,
    .context = NULL
  };

  int N = atoi(argv[2]);
  int i0 = N/2;
  ivec2 shape = {N, N};
  dvec2 xymin = {-1, -1};
  dbl h = 2.0/(N-1);

  /**
   * For the 'v' slowness function, set the domain to be [0, 1] x [0,
   * 1] with a point source at (0, 0).
   */
  if (slo_fun == 'v') {
    i0 = 0;
    xymin = (dvec2) {0, 0};
    h = 1.0/(N-1);
  }

  /**
   * For the 'g' slowness function, set the domain to [0, 1/2] x [0,
   * 1/2] with a point source at (0, 0).
   */
  if (slo_fun == 'g') {
    i0 = 0;
    xymin = (dvec2) {0, 0};
    h = 0.5/(N-1);
  }

  eik_init(scheme, &slow, shape, xymin, h);

  /**
   * Ensure that the factor radius is at least 0.1 in each domain.
   */
  int R = N/20;
  if (slo_fun == 'v') {
    R = N/10;
  }
  if (slo_fun == 'g') {
    R = N/5;
  }
  if (R < 5) R = 5;

  /**
   * Initialize inside box.
   */
  for (int i = 0; i < N; ++i) {
    int abs_di = abs(i - i0);
    dbl x = h*i + xymin.x;
    for (int j = 0; j < N; ++j) {
      int abs_dj = abs(j - i0);
      dbl y = h*j + xymin.y;
      int r = MAX(abs_di, abs_dj);
      if (r <= R) {
        ivec2 ind = {i, j};
        jet_s J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
        if (r < R) {
          eik_add_valid(scheme, ind, J);
        } else {
          eik_add_trial(scheme, ind, J);
        }
      }
    }
  }

  eik_build_cells(scheme);

  eik_solve(scheme);

  jet_s *jets = eik_get_jets_ptr(scheme);

  char path[64];

  FILE *fp = NULL;

  snprintf(path, sizeof(path), "T_%d.bin", N);
  fp = fopen(path, "wb");
  for (int i = 0; i < N*N; ++i) {
    fwrite(&jets[i].f, sizeof(dbl), 1, fp);
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Tx_%d.bin", N);
  fp = fopen(path, "wb");
  for (int i = 0; i < N*N; ++i) {
    fwrite(&jets[i].fx, sizeof(dbl), 1, fp);
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Ty_%d.bin", N);
  fp = fopen(path, "wb");
  for (int i = 0; i < N*N; ++i) {
    fwrite(&jets[i].fy, sizeof(dbl), 1, fp);
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Txy_%d.bin", N);
  fp = fopen(path, "wb");
  for (int i = 0; i < N*N; ++i) {
    fwrite(&jets[i].fxy, sizeof(dbl), 1, fp);
  }
  fclose(fp);

  int k;

  snprintf(path, sizeof(path), "E_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl e = u(x, y) - jets[k++].f;
      fwrite(&e, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Ex_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl ex = ux(x, y) - jets[k++].fx;
      fwrite(&ex, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Ey_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl ey = uy(x, y) - jets[k++].fy;
      fwrite(&ey, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Exx_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl Txx = eik_Txx(scheme, (dvec2) {x, y});
      dbl exx = uxx(x, y) - Txx;
      fwrite(&exx, sizeof(dbl), 1, fp);
      ++k;
    }
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Exy_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl Txy = eik_Txy(scheme, (dvec2) {x, y});
      dbl exy = uxy(x, y) - Txy;
      fwrite(&exy, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

  snprintf(path, sizeof(path), "Eyy_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl Tyy = eik_Tyy(scheme, (dvec2) {x, y});
      dbl eyy = uyy(x, y) - Tyy;
      fwrite(&eyy, sizeof(dbl), 1, fp);
      ++k;
    }
  }
  fclose(fp);

  eik_deinit(scheme);
  eik_dealloc(&scheme);
}
