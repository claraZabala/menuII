#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOL 1E-8

// avalua la funció H en el punt x
void H(const double x[3], double Hx[2]);

// avalua la matriu DH en el punt x.
void DH(const double x[3], double DHx[2][3]);

/*
 * x[] conté un punt conegut de la corba.
 * Si es pot fer correctament la predicció, aquesta es posa a x[] i es retorna
 * el valor 0; si no es pot fer, es retorna el valor 1
 */
int prediccio(int sig, double h, double x[3]);

/*
 * x0[] conté un punt conegut de la corba.
 * Si la correcció funciona bé, el nou punt de la corba és a x[] i es retorna el
 * valor 0; però si Newton no convergeix, es retorna el valor 1
 */
int correccio(double h, double x0[3], double x[3], int kmax, double prec);

/* 
 * Resol sistemes lineals Ax = b de 3 equacions i incògnites. Aquesta funció
 * s’invocarà des de la funció de correcció. Si el determinant de la matriu és
 * no nul (amb una certa tolerància), la solució ha de ser a x[] i cal retornar
 * el valor 0; en cas contrari es retorna el valor 1
 */
int resoldre(double A[3][3], const double b[3], double x[3]);

double norma2(double x[3]);
double determinant(double A[3][3]);

/*
 * (x0)^2 + 1.1(x1)^2 + 0.9(x2)^2 − 0.9 = 0
 * (x0)^2 − 1.2(x1)^2 − (x0 − 1) − x2 = 0
 */
void H(const double x[3], double Hx[2]) {
  Hx[0] = x[0] * x[0] + 1.1 * x[1] * x[1] + 0.9 * x[2] * x[2] - 0.9;
  Hx[1] = x[0] * x[0] - 1.2 * x[1] * x[1] - (x[0] - 1.0) - x[2];
}

/* 
 * Fila 0: derivades d'H respecte x0, x1, x2
 * 2*x0, 2*1.1*x1, 2*0.9*x2
 *
 * Fila 1: derivades d'H respecte x0, x1, x2
 * 2*x0 - 1, -2*1.2*x1, -1
 */
void DH(const double x[3], double DHx[2][3]) {
  DHx[0][0] = 2.0 * x[0];
  DHx[0][1] = 2.2 * x[1];
  DHx[0][2] = 1.8 * x[2];
  DHx[1][0] = 2.0 * x[0] - 1.0;
  DHx[1][1] = -2.4 * x[1];
  DHx[1][2] = -1.0;
}

// Resol un sistema 3x3 usant el mètode de Cramer
int resoldre(double A[3][3], const double b[3], double x[3]) {
  double det, cramer[3][3];
  det = determinant(A);
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
  return 0;
}

// Calcula la norma 2 d'un vector de R3
double norma2(double x[3]) {
  return sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

// Calcula el determinant d'una matriu 3x3
double determinant(double A[3][3]) {
  return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
         A[1][0] * A[2][1] * A[0][2] - A[0][2] * A[1][1] * A[2][0] -
         A[0][1] * A[1][0] * A[2][2] - A[1][2] * A[2][1] * A[0][0];
}

/*
 * FUNCIÓ DE PREDICCIÓ
 * Comprova si DH(P) té rang 2 amb una tolerància.
 * x inicialment és un punt conegut P de la corba i es fa una nova predicció
 * Fer el producte vectorial V = ∇H0(P) × ∇H1(P)
 * El punt predit és Q = P + sig · h ·V / ||V||2
 */
int prediccio(int sig, double h, double x[3]) {
  double DHx[2][3], V[3], normaV, cos;
  DH(x, DHx);
  cos =
      (DHx[0][0] * DHx[1][0] + DHx[0][1] * DHx[1][1] + DHx[0][2] * DHx[1][2]) /
      (norma2(DHx[0]) * norma2(DHx[1]));
  printf("pred: cos=%lf\n", cos);
  if (fabs(cos) < TOL) {
    return 1;
  }
  for (int i = 0; i < 3; i++) {
    V[i] = DHx[0][(i + 1) % 3] * DHx[1][(i + 2) % 3] -
           DHx[0][(i + 2) % 3] * DHx[1][(i + 1) % 3];
  }
  normaV = norma2(V);
  printf("pred: Vt=(%+.2le,%+.2le,%+.2le)\n", V[0] / normaV, V[1] / normaV,
         V[2] / normaV);

  for (int i = 0; i < 3; i++) {
    x[i] += sig * h * V[i] / normaV;
  }
  printf("pred: x=(%+.6le,%+.6le,%+.6le)\n", x[0], x[1], x[2]);
  return 0;
}

/*
 * FUNCIÓ DE CORRECCIÓ 
 * Resol el següent sistema amb Newton-Raphson:
 * H0(x0, x1, x2) = 0
 * H1(x0, x1, x2) = 0
 * (||X − P||2)^2 = h^2   -> aquí X és el punt obtingut a predicció
 *
 * x(k+1) = x(k) + z
 * DF(x(k)) * z = −F(x(k))
 *
 * x0 conté el punt P inicial, x conté el punt predit i es va corregint
 */
int correccio(double h, double x0[3], double x[3], int kmax, double prec) {
  /*
   * F conté el sistema -H0, -H1, h^2-||x-x0||2^2
   * DF la derivada del sistema H0, H1, ||x-x0||2^2 - h^2 respecte de x
   */
  int k = 0;
  double DF[3][3], F[3];
  double Hx[2], DHx[2][3], z[3], err = 1.0;
  while (k < kmax && err > prec) {
    H(x, Hx);
    DH(x, DHx);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        DF[i][j] = DHx[i][j];
      }
      F[i] = -Hx[i];
    }
    for (int i = 0; i < 3; i++) {
      DF[2][i] = 2.0 * (x[i] - x0[i]);
    }
    F[2] = pow(h, 2) - pow(x[0] - x0[0], 2) - pow(x[1] - x0[1], 2) -
           pow(x[2] - x0[2], 2);

    if (resoldre(DF, F, z) == 1) {
      printf("No s'ha pogut resoldre el sistema.\n");
      return 1;
    }

    x[0] += z[0];
    x[1] += z[1];
    x[2] += z[2];
    err = norma2(z);
    k++;
  }

  if (k == kmax && err > prec) {
    printf("corr: err=%+.6le\n", err);
    return 1;
  }

  printf("corr: k = %d, err = %+.6le, x=(%+.6le, %+.6le, %+.6le)\n", k, err,
         x[0], x[1], x[2]);
  return 0;
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
  for (int j = 0; j < 3; j++) {
    fgets(str, 30, stdin);
    x_aux = strtod(str, &ptr);
    x[j] = x_aux;
  }
  if (errno == EINVAL) {
    printf("x introduït incorrectament!");
    return 1;
  }
  printf("Has introduit el punt x=(%lf, %lf, %lf)\n", x[0], x[1], x[2]);
  H(x, Hx);
  printf("Hx = (%lf, %lf)\n", Hx[0], Hx[1]);
  if (Hx[0] > prec || Hx[1] > prec) {
    printf("El x introduït no és de la corba\n");
    return -1;
  }
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