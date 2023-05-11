import numpy as np
from numpy.linalg import norm


def H(x):
    return np.array([x[0] ** 2 + 1.1 * x[1] ** 2 + 0.9 * x[2] ** 2 - 0.9,
                     x[0] ** 2 - 1.2 * x[1] ** 2 - x[0] + 1 - x[2]])


def DH(x):
    return np.array([[2 * x[0], 2.2 * x[1], 1.8 * x[2]], [2 * x[0] - 1, -2.4 * x[1], -1]])


def prediccio(sig, h, p):
    DHx = DH(p)
    cos = np.dot(DHx[0], DHx[1]) / (norm(DHx[0]) * norm(DHx[1]))
    print("pred: cos =", cos)
    if cos == 0:
        return np.zeros(3), False
    v = np.cross(DHx[0], DHx[1])
    normaV = norm(v, 2)
    v_normalitzat = v / normaV
    q = p + sig * h * v_normalitzat
    return q, True


def correccio(h, p, x, kmax, prec):
    k = 0
    error = np.inf
    while k < kmax and error > prec:
        DF = np.append(DH(x), [2 * (x - p)], axis=0)
        F = np.append(H(x), [norm(x - p) ** 2 - h ** 2], axis=0)
        menysF = -F
        det = np.linalg.det(DF)
        c1 = np.linalg.det(np.transpose(np.append(np.append([menysF], [DF[:, 1]], axis=0), [DF[:, 2]], axis=0)))
        c2 = np.linalg.det(np.transpose(np.append(np.append([DF[:, 0]], [menysF], axis=0), [DF[:, 2]], axis=0)))
        c3 = np.linalg.det(np.transpose(np.append(np.append([DF[:, 0]], [DF[:, 1]], axis=0), [menysF], axis=0)))
        z = np.array([c1/det, c2/det, c3/det])
        x += z
        error = norm(z)
        k += 1
    if error > prec:
        return np.zeros(3), False
    print("corr: k =", k, "error =", error, "x =", x)
    return x, True


def main():
    i = 0
    h = 0.01
    N = 650
    kmax = 8
    prec = 1e-8
    p = np.array([0, 0, 1])
    print("Punt ", i, ": ", p)
    sig = 1

    while i < N:
        q, valor = prediccio(sig, h, p)
        if valor:
            p, valor = correccio(h, p, q, kmax, prec)
            if valor:
                i += 1
                print("Punt", i, ":", p)
            else:
                return 1
        else:
            return 1


main()
main()