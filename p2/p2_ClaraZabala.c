#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOL 1E-4

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

double norma2(double x[3]);
double determinant(double A[3][3]);

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
  double det, cramer[3][3];
  det = determinant(A);
  // printf("det=%+.6le\n", det);
  if (fabs(det) < TOL) {
    return 1;
  }
  for (int i = 0; i < 3; i++) {
    cramer[i][0] = b[i];
    for (int j = 1; j < 3; j++) {
      cramer[i][j] = A[i][j];
    }
  }
  x[0] = determinant(cramer) / det;
  for (int i = 0; i < 3; i++) {
    cramer[i][1] = b[i];
    for (int j = 0; j < 3; j += 2) {
      cramer[i][j] = A[i][j];
    }
  }
  x[1] = determinant(cramer) / det;
  for (int i = 0; i < 3; i++) {
    cramer[i][2] = b[i];
    for (int j = 0; j < 2; j++) {
      cramer[i][j] = A[i][j];
    }
  }
  x[2] = determinant(cramer) / det;
  // printf("z = (%+.6le, %+.6le, %+.6le)\n", x[0], x[1], x[2]);
  return 0;
  /* GAUSS - no funciona perquè no comprova si són 0 les components
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 3; j++) {
      c = A[j][i] / A[i][i];
      for (int k = i; k < 4; k++) {
        if (k == 3) {
          b[j] -= c * b[i];
          printf("b[%d] = %f\n", j, b[j]);
        } else {
          A[j][k] -= c * A[i][k];
          printf("A[%d][%d] = %f\n", j, k, A[j][k]);
        }
      }
    }
  }
  x[2] = b[2] / A[2][2];
  for (int i = 1; i > -1; i--) {
    s = 0;
    for (int j = i + 1; j < 3; j++) {
      s += (A[i][j] * x[j]);
      x[i] = (b[i] - s) / A[i][i];
    }
  }
  */
}

double determinant(double A[3][3]) {
  return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
         A[1][0] * A[2][1] * A[0][2] - A[0][2] * A[1][1] * A[2][0] -
         A[0][1] * A[1][0] * A[2][2] - A[1][2] * A[2][1] * A[0][0];
}

/*
 * Comprovar si DH(P) té rang 2 amb una tolerància.
 * x inicialment és un punt conegut P de la corba i es fa una nova predicció
 * Fer el producte vectorial V = ∇H0(P) × ∇H1(P)
 * El punt predit és Q = P + sig · h ·V / ||V||2
 */
int prediccio(int sig, double h, double x[3]) {
  double DHx[2][3], V[3], normaV, cos = 0.0;
  DH(x, DHx);
  for (int i = 0; i < 3; i++) {
    cos += DHx[0][i] * DHx[1][i];
  }
  cos /= (norma2(DHx[0]) * norma2(DHx[1]));
  printf("pred: cos=%lf\n", cos);
  if (fabs(cos) > TOL) {
    V[0] = DHx[0][1] * DHx[1][2] - DHx[0][2] * DHx[1][1];
    V[1] = DHx[0][2] * DHx[1][0] - DHx[0][0] * DHx[1][2];
    V[2] = DHx[0][0] * DHx[1][1] - DHx[0][1] * DHx[1][0];
    normaV = norma2(V);
    printf("pred: Vt=(%+.2le,%+.2le,%+.2le)\n", V[0] / normaV, V[1] / normaV,
           V[2] / normaV);

    for (int i = 0; i < 3; i++) {
      x[i] += sig * h * V[i] / normaV;
    }
    printf("pred: x=(%+.6le,%+.6le,%+.6le)\n", x[0], x[1], x[2]);
    return 0;
  }
  return 1;
}

/* resoldre el sistema amb Newton-Raphson
 H0(x0, x1, x2) = 0
 H1(x0, x1, x2) = 0
 (||X − P||2)^2 = h^2   -> aquí X és el punt obtingut a predicció
  -
  x(k+1) = x(k) + z
  DF(x(k)) * z = −F(x(k))
 */
// x0 conté el punt P inicial, x conté el punt predit i es va corregint
int correccio(double h, double x0[3], double x[3], int kmax, double prec) {
  int k = 0;
  double DF[3][3], F[3];
  double Hx[2], DHx[2][3], z[3], err = 1.0;
  H(x, Hx);
  DH(x, DHx);
  while (k < kmax && err > prec) {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        DF[i][j] = DHx[i][j];
      }
      F[i] = -1.0 * Hx[i];
    }
    for (int i = 0; i < 3; i++) {
      DF[2][i] = 2.0 * (x[i] - x0[i]);
    }
    F[2] = -1.0 * (pow(x[0] - x0[0], 2) + pow(x[1] - x0[1], 2) +
                   pow(x[2] - x0[2], 2) - pow(h, 2));
    /*printf("DF=(%+.2le,%+.2le,%+.2le)\n(%+.2le,%+.2le,%+.2le)\n"
           "(%+.2le,%+.2le,%+.2le),\n-F=(%+.2le,%+.2le,%+.2le)\n",
           DF[0][0], DF[0][1], DF[0][2], DF[1][0], DF[1][1], DF[1][2], DF[2][0],
           DF[2][1], DF[2][2], F[0], F[1], F[2]);
           */
    if (resoldre(DF, F, z) == 1) {
      printf("No s'ha pogut resoldre el sistema.\n");
      return 1;
    }
    x[0] += z[0];
    x[1] += z[1];
    x[2] += z[2];
    err = norma2(z);
    H(x, Hx);
    DH(x, DHx);
    k++;
  }
  if (k == kmax && err > prec) {
    printf("Newton no ha convergit en les %d iteracions.\n", kmax);
    return 1;
  }
  printf("corr: k = %d, err = %+.6le, x=(%+.6le, %+.6le, %+.6le)\n", k, err,
         x[0], x[1], x[2]);
  return 0;
}

double norma2(double x[3]) {
  return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

int main(void) {
  /*
   * h: distància que es vol que hi hagi entre punts consecutius de la corba,
   * N: quantitat de punts de la corba que es volen calcular,
   * prec: precisió que es demana al mètode de Newton entre iterats consecutius,
   * kmax: nombre màxim d’iteracions permeses al mètode de Newton,
   */
  double h = 0.01, prec = 1e-8, x[3], x0[3], x_aux, Hx[2];
  int N = 650, kmax = 8, i = 0;
  long sig;
  char *ptr = NULL, str[30];
  FILE *fitxer;
  fitxer = fopen("punts_corba.txt", "w");
  if (fitxer == NULL) {
    printf("Error creant el fitxer punts_corba.txt");
    return 1;
  }
  /*
   * Llegim un punt inicial x ∈ R que sigui de la corba i un valor sig ∈ {+1,
   * −1}.
   */
  printf("Introdueix un punt inicial x de R3 que sigui de la corba (número + "
         "intro + número + intro + número + intro):\n");
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  x[0] = x_aux;
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  x[1] = x_aux;
  fgets(str, 30, stdin);
  x_aux = strtod(str, &ptr);
  x[2] = x_aux;
  if (errno == EINVAL) {
    printf("x introduït incorrectament!");
    return 1;
  }
  H(x, Hx);
  printf("Hx = (%lf, %lf)\n", Hx[0], Hx[1]);
  if (Hx[0] > prec || Hx[1] > prec) {
    printf("El x introduït no és de la corba\n");
    return -1;
  }
  printf("Has introduit el punt (%lf, %lf, %lf)\n", x[0], x[1], x[2]);
  fprintf(fitxer, "%d %+.6le %+.6le %+.6le\n", i, x[0], x[1], x[2]);
  printf("Introdueix un valor sig que sigui 1 o -1: ");
  fgets(str, 30, stdin);
  sig = strtol(str, &ptr, 10);
  if (sig != 1 && sig != -1) {
    printf("sig incorrecte!");
    return 1;
  }
  while (i < N) {
    // guardem el punt conegut de la corba per la correcció
    x0[0] = x[0];
    x0[1] = x[1];
    x0[2] = x[2];
    printf("punt i = %d\n", i + 1);
    if (prediccio((int)sig, h, x) == 0) {
      if (correccio(h, x0, x, kmax, prec) == 0) {
        i++;
        fprintf(fitxer, "%d %+.6le %+.6le %+.6le\n", i, x[0], x[1], x[2]);
      } else {
        printf("Newton no ha convergit a la correcció del punt %d\n", i + 1);
        return 1;
      }
    } else {
      printf("No s'ha pogut fer bé la predicció del punt %d\n", i + 1);
      return 1;
    }
  }
  fclose(fitxer);
  return 0;
}