#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define A (-0.8)
#define B 0.5
#define C (-1.2)
#define D (-0.3)

#define TOL 1E-6
#define pre 1E-8
#define llindar 1E-3
#define iter 10
#define n 1000

/* Struct que ens servirà per representar els punts */
struct xy {
    double x, y;// coordenades
    int cp;     // control de pas, guarda el nombre d'iteracions
};

/* Funció de la que volem trobar els 0 */
double func_f(double x, double y) {
    double primer, segon;
    primer = 16. * pow(x, 4) + pow(y, 4) + A * x * pow(y, 2) - 16. * pow(x, 2) * y - 1.;
    segon = pow(x, 2) + pow(y - 1., 2) + B * x * y + C;
    return primer * segon + D;
}

/* Derivada x de la f */
double dfx(double x, double y) {
    double primer, segon, dxprimer, dxsegon;
    primer = 16. * pow(x, 4) + pow(y, 4) + A * x * pow(y, 2) - 16. * pow(x, 2) * y - 1.;
    segon = pow(x, 2) + pow(y - 1., 2) + B * x * y + C;
    dxprimer = 64 * pow(x, 3) + A * pow(y, 2) - 32 * x * y;
    dxsegon = 2 * x + B * y;
    return primer * dxsegon + dxprimer * segon;
}

/* Derivada y de la f */
double dfy(double x, double y) {
    double primer, segon, dyprimer, dysegon;
    primer = 16. * pow(x, 4) + pow(y, 4) + A * x * pow(y, 2) - 16. * pow(x, 2) * y - 1.;
    segon = pow(x, 2) + pow(y - 1., 2) + B * x * y + C;
    dyprimer = 4 * pow(y, 3) + 2 * A * x * y - 16 * pow(x, 2);
    dysegon = 2 * (y - 1) + B * x;
    return primer * dysegon + dyprimer * segon;
}

/* Norma 2 entre 2 vectors en dimensió 2 */
double norma2(double a, double b) {
    return sqrt(pow(a, 2) + pow(b, 2));
}

/* Mètode que busca un bon punt inicial.
   També guarda en un document tots els zeros que troba de la funció.
   Es podria fer que retornés més d'un punt inicial. */
struct xy puntsInicials() {
    int num = 10000, i, j;
    struct xy xy0;
    double minim = 1e12, aux, x, y, pas;
    FILE *fitxer;

    fitxer = fopen("zeros.txt", "w+");
    if (fitxer == NULL) {
        printf("Error obrint el fitxer");
        xy0.x = FP_NAN;
        xy0.y = FP_NAN;
        return xy0;
    }

    x = -1.5;
    y = -1.5;
    pas = 8.0 / num;
    printf("Pas %.4f\n", pas);
    for (i = 0; i < num; i++) {
        for (j = 0; j < num; j++) {
            aux = func_f(x, y);
            //Imprimim al fitxer els punts propers a 0.
            if (fabs(aux) < llindar) {
                fprintf(fitxer, "%f %f %f \n", x, y, aux);
            }
            //Alhora busquem el punt que doni un resultat mínim.
            if (fabs(aux) < fabs(minim)) {
                xy0.x = x;
                xy0.y = y;
                minim = aux;
            }
            x += pas;
        }
        x = -1.5;
        y += pas;
    }
    fclose(fitxer);
    return xy0;
}

/* Mètode per fer la correcció del punt inicial fent el mètode de Newton per 1 variable */
struct xy correccioInicial(struct xy xy0) {
    struct xy correc;
    double dx, dy, t = 0.0;
    double g, gprima, x, y;
    int i = 0;
    x = xy0.x;
    y = xy0.y;

    while (i < iter && (fabs(func_f(x, y)) > pre)) {
        dx = dfx(x, y);
        dy = dfy(x, y);
        g = func_f(x, y);
        gprima = pow(dx, 2) + pow(dy, 2);

        t = t - g / gprima;
        x = xy0.x + t * dx;
        y = xy0.y + t * dy;
        printf("t = %.6le, x = %.6le, y = %.6le \n", t, x, y);

        i++;
    }

    correc.x = x;
    correc.y = y;
    printf("Hem fet %d iteracions per corregir la predicció inicial\n", i);
    return correc;
}

/* Mètode per predir el següent punt seguint la tangent del punt anterior */
struct xy prediccio(double x0, double y0, double h, int direccio) {
    struct xy xy1;
    double dx, dy, norma;
    dx = dfx(x0, y0);
    dy = dfy(x0, y0);
    norma = norma2(dx, dy);
    if (norma < TOL) {
        printf("Som molt a prop d'un punt singular, aturem la continuació");
        xy1.x = FP_NAN;
        xy1.y = FP_NAN;
    } else {
        xy1.x = x0 + direccio * h * dy / norma;
        xy1.y = y0 - direccio * h * dx / norma;
    }
    return xy1;
}

/* Mètode per corretgir un punt usant el mètode de Newton en dues dimensions */
struct xy correccio(struct xy xy0, struct xy xy1, double h) {
    double f2, a, b, c, d, xk, yk, det;
    struct xy xyk1;
    int iteracio = 0;
    xk = xy1.x;
    yk = xy1.y;

    while (iteracio < iter && (fabs(func_f(xk, yk)) > pre)) {
        f2 = pow(xk - xy0.x, 2) + pow(yk - xy0.y, 2) - pow(h, 2);
        a = dfx(xk, yk);
        b = dfy(xk, yk);
        c = 2 * (xk - xy0.x);
        d = 2 * (yk - xy0.y);
        det = a * d - c * b;

        if (fabs(det) > pre) {
            xk = xk - (d * func_f(xk, yk) - b * f2) / det;
            yk = yk - (a * f2 - c * func_f(xk, yk)) / det;
            iteracio++;
        } else {
            printf("!!! Determinant massa proper a 0 a la correció %d\n", iteracio);
            break;
        }
    }

    xyk1.x = xk;
    xyk1.y = yk;
    xyk1.cp = iteracio;

    if (iteracio == iter) {
        printf("No hem convergit en la correcció\n");
    } else {
        //printf("Nombre d'iteracions correcció: %d\n", iteracio);
    }

    return xyk1;
}

/* Mètode que ajunta els dos anteriors i va gestionant les iteracions per trobar tots els punts de la corba */
int predictorCorrector(struct xy puntIni, FILE *fitxer, int direccio) {
    struct xy xy0, xy1;
    int punts = 0, i;
    double h = 1E-2;
    int comptador[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    //Al comptador comptem quantes vegades s'ha tardat x nombre d'iteracions de Newton per arribar al resultat

    printf("Comencem el mètode pel punt (%f, %f) en la direcció %d\n", puntIni.x, puntIni.y, direccio);

    xy0.x = puntIni.x;
    xy0.y = puntIni.y;

    fprintf(fitxer, "%f %f\n", xy0.x, xy0.y);

    while (punts < n) {
        xy1 = prediccio(xy0.x, xy0.y, h, direccio);
        if (isnan(xy1.x) || isnan(xy1.y)) {
            printf("Aturem el programa\n");
            return -1;
        } else {
            xy0 = correccio(xy0, xy1, h);
            comptador[xy0.cp] += 1;
            if (xy0.cp > 3) {
                h = h / 2;
            } else if (xy0.cp < 3) {
                h = h * 1.1;
            }
            if (fabs(func_f(xy0.x, xy0.y)) < 0.2) {
                fprintf(fitxer, "%f %f\n", xy0.x, xy0.y);
                punts++;
            } else {
                printf("Aturem perquè ens hem anat de la corba\n");
                return -1;
            }
        }
    }
    fprintf(fitxer, "\n");

    printf("Nombre d'iteracions fetes per cada correció:\n");
    printf("[");
    for (i = 0; i < 10; i++) {
        printf("%d:%d ", i, comptador[i]);
    }
    printf("]\n");
    return 0;
}

/* Mètode que, a partir d'uns punts inicials, aplicarà el corrector predictor en ambdues direccions*/
void calcularPunts(struct xy punts[], int length, FILE *fitxer) {
    struct xy puntIni;
    int i;

    for (i = 0; i < length; i++) {
        puntIni = punts[i];
        if (predictorCorrector(puntIni, fitxer, 1) == 0) {
            printf("Mètode completat correctament.\n");
        }

        if (predictorCorrector(puntIni, fitxer, -1) == 0) {
            printf("Mètode completat correctament.\n");
        }
    }
}

int main(void) {
    struct xy puntIni, punts[2];// punts permet fins 2 punts inicials
    long metode;
    char *ptr, str[30];
    FILE *fitxer;

    // Possibilitat de demanar els paràmetres per pantalla (A,B entre -1 i 1, C,D entre -2 i 0)
    printf("Mètode de continuació, paràmetres A = %.1f, B = %.1f, C = %.1f, D = %.1f \n", A, B, C, D);
    fitxer = fopen("continuacio.txt", "w+");
    if (fitxer == NULL) {
        printf("Error obrint el fitxer");
        return -1;
    }

    printf("Calculem el punt inicial\n");
    printf("Vols usar un punt inicial precalculat? [0=Sí, 1=No]: ");
    fgets(str, 30, stdin);
    metode = strtol(str, &ptr, 10);

    if (metode == 0) {
        //Usarem punts inicials ja calculats
        // A = -0.8, B = 0.5, C = -1.2, D = -0.3;
        // 0.806000 1.162000 punt inicial de dins
        // 0.517600 1.803900 punt inicial de fora

        // A = -0.8, B = 0.5, C = -1.2, D = 0.3;
        //x = -0.198000, y = -0.726000 Component de baix
        //x = -0.964000, y = 2.104000 Component de l'esquerra
        //x = 0.704000, y = 1.488000 Component de la dreta

        punts[0].x = 0.806000;
        punts[0].y = 1.162000;
        punts[1].x = 0.517600;
        punts[1].y = 1.803900;

        calcularPunts(punts, 2, fitxer);
    } else {
        // Buscarem un punt inicial adient
        puntIni = puntsInicials();
        printf("El bon punt inicial trobat és x = %.6f, y = %.6f \n", puntIni.x, puntIni.y);
        printf("Fem la correcció del punt inicial\n");
        puntIni = correccioInicial(puntIni);
        printf("El punt inicial corregit és x = %.6f, y = %.6f \n", puntIni.x, puntIni.y);

        punts[0] = puntIni;
        calcularPunts(punts, 1, fitxer);
    }

    fclose(fitxer);
    printf("S'ha acabat el mètode, els punts es troben a continuacio.txt\n");

    return 0;
}
