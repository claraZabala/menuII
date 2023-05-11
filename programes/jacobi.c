#include<stdio.h>
#include<math.h>

/* Arrange systems of linear
   equations to be solved in
   diagonally dominant form
   and form equation for each
   unknown and define here
*/
/* In this example we are solving
   4x - y - z = 1
   -x + 4y - t = 2
   -x + 4z - t = 0
	 - y - z + 4t = 1
*/
/* Equations:
   x = (1+y+z)/4
   y = (2+x+t)/4
   z = (x+t)/4
	 t = (1+y+z)/4
*/
/* Defining function */
#define f1(x,y,z,t)  (1+y+z)/4.0
#define f2(x,y,z,t)  (2+x+t)/4.0
#define f3(x,y,z,t)  (x+t)/4.0
#define f4(x,y,z,t)	 (1+y+z)/4.0

/* Main function */
int main()
{
 float x0=0.0, y0=0.0, z0=0.0, t0=0.0, x1, y1, z1, t1, e1, e2, e3, e4, e;
 int count=1;
 printf("Entra l'error de toleràcia:\n");
 scanf("%f", &e);

 printf("\nIter\tx\ty\tz\tt\n");
 do
 {
  /* Calculation */
  x1 = f1(x0,y0,z0,t0);
  y1 = f2(x0,y0,z0,t0);
  z1 = f3(x0,y0,z0,t0);
	t1 = f4(x0,y0,z0,t0);
  printf("%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n", count, x1,y1,z1,t1);

  /* Error */
  e1 = fabs(x0-x1);
  e2 = fabs(y0-y1);
  e3 = fabs(z0-z1);
	e4 = fabs(t0-t1);

  count++;

  /* Set value for next iteration */
  x0 = x1;
  y0 = y1;
  z0 = z1;
	t0 = t1;

 }while(e1>e || e2>e || e3>e || e4>e);

 printf("\nSolució: x=%0.3f, y=%0.3f, z=%0.3f i t=%0.3f\n",x1,y1,z1,t1);
 return 0;
}