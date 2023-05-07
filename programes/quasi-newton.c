#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Defining equation to be solved.
   Change this equation to solve another problem. */
int F(double x, double y, double F[2]) {
  F[0] = x * x * x + y * y + x + y - 2.0;
  F[1] = x * y * y + x * x + x - y;
}

/* Defining derivative of f(x).
   As you change f(x), change this function also. */
int DF(double x, double y, double DF[2][2]) {
  DF[0][0] = 3.0 * x * x + 1.0;
  DF[0][1] = 2.0 * y + 1.0;
  DF[1][0] = y * y + 2.0 * x + 1.0;
  DF[1][1] = 2.0 * x * y - 1.0;
  return 0;
}

int main() {
  float e;
  double DFx[2][2], Fx[2], x0, x1, y0, y1, det, err;
  int step = 1, N;
  printf("\nEnter initial guess:\n");
  scanf("%lf %lf", &x0, &y0);
  printf("Enter tolerable error:\n");
  scanf("%f", &e);
  printf("Enter maximum iteration:\n");
  scanf("%d", &N);

  DF(x0, y0, DFx);

  do {
    printf("DF=((%lf,%lf)\n(%lf,%lf))\n", DFx[0][0], DFx[0][1], DFx[1][0],
           DFx[1][1]);
    F(x0, y0, Fx);
    det = DFx[0][0] * DFx[1][1] - DFx[0][1] * DFx[1][0];
    // printf("det=%lf\n", det);
    if (fabs(det) < e) {
      printf("Mathematical Error.");
      exit(1);
    }

    x1 = x0 - (DFx[1][1] * Fx[0] - DFx[0][1] * Fx[1]) / det;
    y1 = y0 - (-DFx[1][0] * Fx[0] + DFx[0][0] * Fx[1]) / det;

    printf("%d\t\t%f\t%f\t%f\t%f\n", step, x0, y0, x1, y1);
    err = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    printf("err=%+.2le\n", err);
    x0 = x1;
    y0 = y1;

    step++;

    if (step > N) {
      printf("Not Convergent.");
      exit(1);
    }
  } while (err > e);

  printf("\nRoot is: %lf, %lf", x1, y1);
  return 0;
}
