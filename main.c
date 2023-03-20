#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define p1 3.0
#define q1 2.0
#define r1 1.0
#define p2(x) (1.0 / (x))
#define q2 0.0
#define r2(x) (1.0 / ((x) * (x)))

void troba_vaps(int exemple, int n, double h, double *a, double *b, double *c);

void compute_a(int exemple, int n, double h, int a, double *as);
void compute_b(int exemple, int n, double h, int a, double *b);
void compute_c(int exemple, int n, double h, int a, double *c);

double cota_inf(int n, const double *a, double *b, double *c);
double cota_sup(int n, const double *a, double *b, double *c);

double compute_det(const double *a, const double *b, const double *c,
                   double lambda, int n);
void factoritzacioLU(int n, double t, const double *a, const double *b,
                     const double *c, double *l, double *d);

double potencia_id(int n, double *l, double *d, double *b);
void calcul_y_z(int n, const double *l, const double *d, const double *b,
                const double *x, double *z);
double max_abs_llista(int n, double *a);

void precisio(int exemple, double lambda, int n);
/*
 * as = h * h * (2.0 / (double)(h * h) - 2.0);
 */
void compute_a(int exemple, int n, double h, int a, double *as) {
  double x, ex1;
  ex1 = 2.0 / r1 - pow(h, 2) * q1 / r1;
  for (int i = 0; i < n; i++) {
    switch (exemple) {
    case 1:
      as[i] = ex1;
      break;
    case 2:
      x = (double)a + i * h;
      as[i] = 2.0 / r2(x) - pow(h, 2) * q2 / r2(x);
      break;
    default:
      printf("Nombre d'exemple incorrecte.\n");
    }
  }
}
/*
 * bs = h * h * (-1.0 / (double)(h * h) + 3.0 / (2.0 * h));
 */
void compute_b(int exemple, int n, double h, int a, double *bs) {
  double x, ex1;
  ex1 = -1.0 / r1 - (h * p1) / (2.0 * r1);
  for (int i = 0; i < n - 1; i++) {
    switch (exemple) {
    case 1:
      bs[i] = ex1;
      break;
    case 2:
      x = a + i * h;
      bs[i] = -1.0 / r2(x) - (h * p2(x)) / (2.0 * r2(x));
      break;
    default:
      printf("Nombre d'exemple incorrecte.\n");
    }
  }
  bs[n - 1] = 0.0;
}
/*
 * cs = h * h * (-1.0 / (double)(h * h) - 3.0 / (2.0 * h));
 */
void compute_c(int exemple, int n, double h, int a, double *c) {
  double x, ex1;
  ex1 = -1.0 / r1 + (h * p1) / (2.0 * r1);
  c[0] = 0.0;
  for (int i = 1; i < n; i++) {
    switch (exemple) {
    case 1:
      c[i] = ex1;
      break;
    case 2:
      x = a + i * h;
      c[i] = -1.0 / (double)r2(x) + (h * p2(x)) / (2.0 * (double)r2(x));
      break;
    default:
      printf("Nombre d'exemple incorrecte.\n");
    }
  }
}

double cota_inf(int n, const double *a, double *b, double *c) {
  double m, cota;
  m = a[0] - fabs(b[0]) - fabs(c[0]);
  for (int i = 1; i < n; i++) {
    cota = a[i] - fabs(b[i]) - fabs(c[i]);
    if (cota < m) {
      m = cota;
    }
  }
  return m;
}

double cota_sup(int n, const double *a, double *b, double *c) {
  double m, cota;
  m = a[0] + fabs(b[0]) + fabs(c[0]);
  for (int i = 1; i < n; i++) {
    cota = a[i] + fabs(b[i]) + fabs(c[i]);
    if (cota > m) {
      m = cota;
    }
  }
  return m;
}

double compute_det(const double *a, const double *b, const double *c,
                   double lambda, int n) {
  /*
   * Calcular determinants(d's):
   * d0=1,
   * d1=a1-lambda,
   * di=(ai-lambda)(di-1)-(bi-1)ci(di-2)) -> 3 doubles dn=det(A-lambdaI)
   */
  double d2 = 1.0, d1, d;
  d1 = a[0] - lambda;
  for (int i = 1; i < n; i++) {
    d = (a[i] - lambda) * d1 - b[i - 1] * c[i] * d2;
    d2 = d1;
    d1 = d;
  }
  // printf("Det(%.6le): %.6le\n",lambda/pow(1.0/11.0,2), d);
  return d;
}

void factoritzacioLU(int n, double t, const double *a, const double *b,
                     const double *c, double *l, double *d) {
  d[0] = a[0] - t;
  // printf("lambda=%.6le, a[0]=%.6le\n",lambda,a[0]);
  for (int i = 1; i < n; i++) {
    l[i] = c[i] / d[i - 1];
    d[i] = (a[i] - t) - l[i] * b[i - 1];
  }
}

void calcul_y_z(int n, const double *l, const double *d, const double *b,
                const double *x, double *z) {
  double *y = (double *)malloc(sizeof(double) * n);
  // Calcul vector y: L*y=xk-1
  y[0] = x[0];
  for (int i = 1; i < n; i++) {
    y[i] = x[i] - l[i] * y[i - 1];
  }
  // Calcul vector z: U*zk=y
  z[n - 1] = y[n - 1] / d[n - 1];
  for (int i = n - 2; i > -1; i--) {
    z[i] = (y[i] - b[i] * z[i + 1]) / d[i];
  }
  free(y);
}

double max_abs_llista(int n, double *a) {
  double max = a[0];
  for (int i = 1; i < n; i++) {
    if (fabs(max) < fabs(a[i])) {
      max = fabs(a[i]);
    }
  }
  return max;
}

double potencia_id(int n, double *l, double *d, double *b) {
  bool no_convergim = true;
  int max_iter = 50, iter = 0, pos;
  double *x = (double *)malloc(sizeof(double) * n);
  double *z = (double *)malloc(sizeof(double) * n);
  double s, s_nova;

  // Creem x0 amb valors entre [-1,1] i un és 1
  // pos = rand() % n;
  pos = 0;
  x[pos] = 1.0;
  for (int i = 0; i < n; i++) {
    if (i != pos) {
      x[i] = (double)random() / RAND_MAX * 2.0 - 1.0;
      // printf("x[%d]=%.6le\n", i, x[i]);
    }
  }
  /*
  bucle per trobar la s correcta (fins que una s i l'anterior siguin molt
  pròximes
  */
  while (no_convergim && iter < max_iter) {
    // Càlcul L*y=xk-1 i U*zk=y
    calcul_y_z(n, l, d, b, x, z);
    s_nova = max_abs_llista(n, z);
    // printf("s_nova: %.6le\n", s_nova);
    if (fabs(s_nova - s) > 1e-3) {
      for (int i = 0; i < n; i++) {
        x[i] = z[i] / s_nova;
      }
    } else {
      no_convergim = false;
    }
    s = s_nova;
    iter++;
  }
  if (iter == max_iter) {
    printf("No hem convergit :(\n");
  }
  free(x);
  free(z);
  return s_nova;
}

void precisio(int exemple, double lambda, int n) {
  double vap;
  switch (exemple) {
  case 1:
    vap = (1.0 + 4.0 * pow(n * M_PI, 2)) / 4.0;
    printf("precisió: %le\n", fabs(vap - lambda));
    break;
  case 2:
    vap = pow((n * M_PI / log(2)), 2);
    printf("precisió: %le\n", fabs(vap - lambda));
    break;
  default:
    printf("Exemple incorrecte.\n");
  }
}

void troba_vaps(int exemple, int n, double h, double *a, double *b, double *c) {
  double det = 0.0, nou_det;
  double s, lambda, lambda_final;
  int j = 0, count_vap = 1;
  double sm, sM, mu;
  double *l = (double *)malloc(sizeof(double) * n);
  double *d = (double *)malloc(sizeof(double) * n);

  printf("Comencem a calcular l'exemple %d\n", exemple);

  // Calcular sm = min(ai - |bi| - |ci|)
  // Calcular sM = max(ai + |bi| + |ci|)
  if (exemple == 1) {
    sm = a[1] - fabs(b[1]) - fabs(c[1]);
    sM = a[1] + fabs(b[1]) + fabs(c[1]);
  } else {
    sm = cota_inf(n, a, b, c);
    sM = cota_sup(n, a, b, c);
  }
  printf("Interval on hi ha els vaps: [%.6le, %.6le]\n", sm, sM);

  /*
  Fixar s>0 petit i avaluar det(A-lambdaI) per lambda=sm+js, j>=0, fins
  arribar a sM. Quan es detecti un canvi de signe entre un càlcul del det i el
  següent tenim un VAP lambda -> usar mètode potència inv.despl.
  */
  s = 0.001;
  lambda = sm;
  while (lambda <= sM) {
    lambda = sm + j * s;
    nou_det = compute_det(a, b, c, lambda, n);
    if (nou_det * det < 0) {
      printf("Canvi de signe entre:\n");
      printf("s = %.6le \t det(s) = %.6le\ns = %.6le \t det(s) = %.6le\n",
             sm + (j - 1) * s, det, lambda, nou_det);
      printf("%d. VAP trobat: %.6le\n", count_vap, lambda / (h * h));
      factoritzacioLU(n, lambda, a, b, c, l, d);
      mu = potencia_id(n, l, d, b);
      printf("Convergència de pot. inv. desp. a mu = %.8le\n", mu);
      lambda_final = (lambda + 1.0 / mu) / (h * h);
      printf("VAP %d: %.8le\n", count_vap, lambda_final);
      precisio(exemple, lambda_final, count_vap);
      count_vap++;
    }
    det = nou_det;
    j++;
  }

  free(d);
  free(l);
}

int main() {
  int a, b;
  double h;
  char *ptr;
  long n, exemple;
  char str[30];
  double *as, *bs, *cs;

  printf("Exemple 1 o 2?: ");
  fgets(str, 30, stdin);
  exemple = strtol(str, &ptr, 10);
  printf("Introdueix un enter n > 1: ");
  fgets(str, 30, stdin);
  n = strtol(str, &ptr, 10);
  if (n > 1) {
    printf("OK\n");
  } else {
    printf("No és un enter o no és > 1.\n");
    return 0;
  }
  as = (double *)malloc(sizeof(double) * n);
  bs = (double *)malloc(sizeof(double) * n);
  cs = (double *)malloc(sizeof(double) * n);
  // Definim la seed
  struct timespec ts;
  srandom(ts.tv_nsec ^ ts.tv_sec); /* Seed the PRNG */
  switch (exemple) {
  case 1:
    a = 0;
    b = 1;
    h = (double)(b - a) / (double)(n + 1);
    // Calcular a's (vector 1,...,n)
    // Calcular b's (vector 1,...,n-1) bn=0
    // Calcular c's (vector 2,...,n -> n-1) c1=0
    compute_a(1, (int)n, h, a, as);
    compute_b(1, (int)n, h, a, bs);
    compute_c(1, (int)n, h, a, cs);
    printf("ai: %.4le, bi: %.4le, ci: %.4le\n", as[1], bs[1], cs[1]);
    troba_vaps(1, (int)n, h, as, bs, cs);
    break;
  case 2:
    a = 1;
    b = 2;
    h = (double)(b - a) / (double)(n + 1);
    // Calcular a's (vector 1,...,n)
    // Calcular b's (vector 1,...,n-1) bn=0
    // Calcular c's (vector 2,...,n -> n-1) c1=0
    compute_a(2, (int)n, h, a, as);
    compute_b(2, (int)n, h, a, bs);
    compute_c(2, (int)n, h, a, cs);
    for (int i = 0; i < n; i++) {
      printf("a[%d]=%.4le, b[%d]=%.4le, c[%d]=%.4le\n", i, as[i], i, bs[i], i,
             cs[i]);
    }
    troba_vaps(2, (int)n, h, as, bs, cs);
    break;
  default:
    printf("Nombre de l'exemple erroni\n");
  }
  free(as);
  free(bs);
  free(cs);
  return 0;
}
