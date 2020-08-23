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
  dbl (*uxy)(dbl, dbl);

  if (slo_fun == '1') {
    s = s_1;
    grad_s = grad_s_1;
    u = u_1;
    ux = ux_1;
    uy = uy_1;
    uxy = uxy_1;
  } else if (slo_fun == 'p') {
    s = s_p;
    grad_s = grad_s_p;
    u = u_p;
    ux = ux_p;
    uy = uy_p;
    uxy = uxy_p;
  } else if (slo_fun == 'v') {
    s = s_v;
    grad_s = grad_s_v;
    u = u_v;
    ux = ux_v;
    uy = uy_v;
    uxy = uxy_v;
  // } else if (slo_fun == 'm') {
  //   s = s_m;
  //   grad_s = grad_s_m;
  //   u = u_m;
  //   ux = ux_m;
  //   uy = uy_m;
  //   uxy = uxy_m;
  // } else if (slo_fun == 'g') {
  //   s = s_g;
  //   grad_s = grad_s_g;
  //   u = u_g;
  //   ux = ux_g;
  //   uy = uy_g;
  //   uxy = uxy_g;
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
  eik_init(scheme, &slow, shape, xymin, h);

  int R = N/20;
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

  /**
   * Initialize inside disk.
   */
  // for (int i = 0; i < N; ++i) {
  //   int di = i - i0, di_sq = di*di;
  //   for (int j = 0; j < N; ++j) {
  //     int dj = j - i0, dj_sq = dj*dj;
  //     dbl r = sqrt(di_sq + dj_sq);
  //     if (r < R) {
  //       dbl x = h*i + xymin.x;
  //       dbl y = h*j + xymin.y;
  //       jet J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
  //       eik_add_valid(scheme, (ivec2) {i, j}, J);
  //     }
  //   }
  // }
  // {
  //   dbl di[4] = {1, 0, -1,  0};
  //   dbl dj[4] = {0, 1,  0, -1};
  //   for (int i = 0; i < N; ++i) {
  //     for (int j = 0; j < N; ++j) {
  //       ivec2 ind = {i, j};
  //       if (eik_get_state(scheme, ind) != VALID) {
  //         continue;
  //       }
  //       for (int k = 0; k < 4; ++k) {
  //         int i_ = i + di[k];
  //         int j_ = j + dj[k];
  //         if (0 <= i_ && i_ < N && 0 <= j_ && j_ < N) {
  //           ivec2 ind_ = {i_, j_};
  //           state_e state = eik_get_state(scheme, ind_);
  //           if (state != VALID && state != TRIAL) {
  //             dbl x = h*i_ + xymin.x;
  //             dbl y = h*j_ + xymin.y;
  //             jet J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
  //             eik_add_trial(scheme, ind_, J);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

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

  snprintf(path, sizeof(path), "Exy_%d.bin", N);
  fp = fopen(path, "wb");
  k = 0;
  for (int i = 0; i < N; ++i) {
    dbl x = xymin.x + h*i;
    for (int j = 0; j < N; ++j) {
      dbl y = xymin.y + h*j;
      dbl exy = uxy(x, y) - jets[k++].fxy;
      fwrite(&exy, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

  eik_deinit(scheme);
  eik_dealloc(&scheme);
}
