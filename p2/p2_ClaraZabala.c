#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOL 1E-8

// avalua la funció H en el punt x
void H(double x[3], double Hx[2]);

// avalua la matriu DH en el punt x.
void DH(double x[3], double DHx[2][3]);

// x[] conté un punt conegut de la corba.
// Si es pot fer correctament la predicció, aquesta es posa a x[] i es retorna
// el valor 0; si no es pot fer, es retorna el valor 1
int prediccio(int sig, double h, double x[3]);

// x0[] conté la predicció.
// Si la correcció funciona bé, el nou punt de la corba és a x[] i es retorna el
// valor 0; però si Newton no convergeix, es retorna el valor 1
int correccio(double h, double x0[3], double x[3], int kmax, double prec);

// Resol sistemes lineals Ax = b de 3 equacions i incògnites. Aquesta funció
// s’invocarà des de la funció de correcció. Si el determinant de la matriu és
// no nul (amb una certa tolerància), la solució ha de ser a x[] i cal retornar
// el valor 0; en cas contrari es retorna el valor 1
int resoldre(double A[3][3], double b[3], double x[3]);

// (x0)^2 + 1.1(x1)^2 + 0.9(x2)^2 − 0.9 = 0
// (x0)^2 − 1.2(x1)^2 − (x0 − 1) − x2 = 0
void H(double x[3], double Hx[2]) {
  Hx[0] = x[0] * x[0] + 1.1 * x[1] * x[1] + 0.9 * x[2] * x[2] - 0.9;
  Hx[1] = x[0] * x[0] - 1.2 * x[1] * x[1] - (x[0] - 1.0) - x[2];
}

// Fila 0: derivades d'H respecte x0, x1, x2
// Fila 1: derivades d'H respecte x0, x1, x2
void DH(double x[3], double DHx[2][3]) {
  DHx[0][0] = 2.0 * x[0];
  DHx[0][1] = 2.0 * 1.1 * x[1];
  DHx[0][2] = 2.0 * 0.9 * x[2];
  DHx[1][0] = 2.0 * x[0] - 1.0;
  DHx[1][1] = 2.0 * 1.2 * x[1];
  DHx[1][2] = -1.0;
}

int resoldre(double A[3][3], double b[3], double x[3]) {
  double det, A_inv[3][3];
  det = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
        A[1][0] * A[2][1] * A[0][2] - A[0][2] * A[1][1] * A[2][0] -
        A[0][1] * A[1][0] * A[2][2] - A[1][2] * A[2][1] * A[0][0];
  if (fabs(det) < TOL) {
    return 1;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A_inv[i][j] =
          ((A[(i + 1) % 3][(j + 1) % 3] * A[(i + 2) % 3][(j + 2) % 3]) -
           (A[(i + 1) % 3][(j + 2) % 3] * A[(i + 2) % 3][(j + 1) % 3])) /
          det;
    }
    x[i] = A[i][0] * b[0] + A[i][1] * b[1] + A[i][2] * b[2];
  }
  return 0;
}

int prediccio(int sig, double h, double x[3]) { return 0; }

/* resoldre el sistema amb Newton-Raphson
 H0(x0, x1, x2) = 0
 H1(x0, x1, x2) = 0
 (||X − P||2)^2 = h^2   -> aquí X és el punt obtingut a predicció
  -
  x(k+1) = x(k) + z
  DF(x(k)) * z = −F(x(k))
 */
// x0 conté el punt inicial, x el punt predit
int correccio(double h, double x0[3], double x[3], int kmax, double prec) {
  int k = 0;
  double DF[3][3], F[3];
  double Hx[2], DHx[2][3], z[3];
  x[0] = x0[0];
  x[1] = x0[1];
  x[2] = x0[2];
  H(x, Hx);
  while (k < kmax && (fabs(Hx[0]) > prec || fabs(Hx[1]) > prec)) {
    DH(x, DHx);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        DF[i][j] = DHx[i][j];
      }
      F[0] = -1.0 * Hx[0];
      F[1] = -1.0 * Hx[1];
    }
    DF[2][0] = 2.0 * (x[0] - x0[0]);
    DF[2][1] = 2.0 * (x[1] - x0[1]);
    DF[2][2] = 2.0 * (x[2] - x0[2]);
    F[2] = -1.0 * (pow(x[0] - x0[0], 2) + pow(x[1] - x0[1], 2) +
                   pow(x[2] - x0[2], 2) - h * h);
    if (resoldre(DF, F, z) == 1) {
      return 1;
    }
    x[0] = x[0] + z[0];
    x[1] = x[1] + z[1];
    x[2] = x[2] + z[2];
    H(x, Hx);
    k++;
  }
  if (k <= kmax) {
    return 0;
  }
  return 1;
}

int main(void) {
  /*
   * h: distància que es vol que hi hagi entre punts consecutius de la corba,
   * N: quantitat de punts de la corba que es volen calcular,
   * prec: precisió que es demana al mètode de Newton entre iterats consecutius,
   * kmax: nombre màxim d’iteracions permeses al mètode de Newton,
   */
  double h = 0.01, prec = 1e-8, x[3], p[3], x_aux;
  int N = 650, kmax = 8, i = 0, seguim = 0;
  long sig;
  char *ptr = NULL, str[30];
  FILE *fitxer;
  fitxer = fopen("punts_corba.txt", "w+");
  if (fitxer == NULL) {
    printf("Error creant el fitxer punts_corba.txt");
    return -1;
  }
  /*
   * Llegim un punt inicial x ∈ R que sigui de la corba i un valor sig ∈ {+1,
   * −1}.
   */
  printf("Introdueix un punt inicial P de R3 que sigui de la corba (número + "
         "intro + número + intro + número + intro):\n");
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  p[0] = x_aux;
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  p[1] = x_aux;
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  p[2] = x_aux;
  if (errno == EINVAL) {
    printf("P introduït incorrectament!");
    return -1;
  }
  printf("Has introduit el punt (%lf, %lf, %lf)\n", p[0], p[1], p[2]);
  fprintf(fitxer, "%d %lf %lf %lf\n", i, p[0], p[1], p[2]);
  printf("Introdueix un valor sig que sigui 1 o -1: ");
  fgets(str, 30, stdin);
  sig = strtol(str, &ptr, 10);
  if (sig == 1 || sig == -1) {
    while (i < N && seguim == 0) {
      if (prediccio((int)sig, h, x) == 0) {
        if (correccio(h, p, x, kmax, prec) == 0) {
          i++;
          fprintf(fitxer, "%d %lf %lf %lf\n", i, x[0], x[1], x[2]);
        } else {
          printf("Newton no ha convergit a la correcció del punt %d\n", i);
          seguim = 1;
        }
      } else {
        printf("No s'ha pogut fer bé la predicció del punt %d\n", i);
        seguim = 1;
      }
    }
  } else {
    printf("sig incorrecte!");
    return -1;
  }
}