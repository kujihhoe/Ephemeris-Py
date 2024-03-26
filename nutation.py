from math import cos, floor, sin
import numpy as np
from numba import jit
from nutation_list import nals_t, cls_t, napl_t, cpl_t


pi2 = 6.283185307179586476925287
S2R = 4.848136811095359935899141e-6
MMS2R = 4.848136811095359935899141e-13
c = 299792.458


@jit
def nutaMx(T):
    def fmod(X, m):
        return X - floor(X / m) * m

    def r1(a):
        return np.array([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]])

    def r3(g):
        return np.array([[cos(g), sin(g), 0], [-sin(g), cos(g), 0], [0, 0, 1]])

    def nutation(T):
        def aa(T):
            F = [
                134.96340251 * 3600
                + 1717915923.2178 * T
                + 31.8792 * T**2
                + 0.051635 * T**3
                - 0.00024470 * T**4,
                357.52910918 * 3600
                + 129596581.0481 * T
                - 0.5532 * T**2
                + 0.000136 * T**3
                - 0.00001149 * T**4,
                93.27209062 * 3600
                + 1739527262.8478 * T
                - 12.7512 * T**2
                - 0.001037 * T**3
                + 0.00000417 * T**4,
                297.85019547 * 3600
                + 1602961601.2090 * T
                - 6.3706 * T**2
                + 0.006593 * T**3
                - 0.00003169 * T**4,
                125.04455501 * 3600
                - 6962890.5431 * T
                + 7.4722 * T**2
                + 0.007702 * T**3
                - 0.00005939 * T**4,
            ]
            F = [fmod(value * S2R, pi2) for value in F]
            return F

        def aa1(T):
            F = [
                2.35555598 + 8328.6914269554 * T,
                6.24006013 + 628.301955 * T,
                1.627905234 + 8433.466158131 * T,
                5.198466741 + 7771.3771468121 * T,
                2.18243920 - 33.757045 * T,
                4.402608842 + 2608.7903141574 * T,
                3.176146697 + 1021.3285546211 * T,
                1.753470314 + 628.3075849991 * T,
                6.203480913 + 334.0612426700 * T,
                0.599546497 + 52.9690962641 * T,
                0.874016757 + 21.3299104960 * T,
                5.481293872 + 7.4781598567 * T,
                5.311886287 + 3.8133035638 * T,
                0.02438175 * T + 0.00000538691 * T**2,
            ]
            for i in range(0, len(F) - 1, 1):
                F[i] = fmod(F[i], pi2)
            return F

        a = aa(T)
        dp = de = dp1 = de1 = 0.0
        for i in range(677, -1, -1):
            arg = sum(nals_t[i][j] * a[j] for j in range(5)) % pi2
            sarg = sin(arg)
            carg = cos(arg)
            dp += (cls_t[i][0] + cls_t[i][1] * T) * sarg + cls_t[i][2] * carg
            de += (cls_t[i][3] + cls_t[i][4] * T) * carg + cls_t[i][5] * sarg
        dpsils = dp
        depsls = de
        a1 = aa1(T)
        for i in range(686, -1, -1):
            arg = sum(napl_t[i][j] * a1[j] for j in range(14)) % pi2
            sarg = sin(arg)
            carg = cos(arg)
            dp1 += cpl_t[i][0] * sarg + cpl_t[i][1] * carg
            de1 += cpl_t[i][2] * sarg + cpl_t[i][3] * carg
        dpsipl = dp1
        depspl = de1
        dpsi = dpsipl + dpsils
        deps = depspl + depsls
        return {"NutaEclp": dpsi * MMS2R, "NutaObliq": deps * MMS2R}

    def obliqAvg(T):
        return (
            84381.406
            - 46.836769 * T
            - 0.0001831 * T**2
            + 0.00200340 * T**3
            - 0.000000576 * T**4
            - 0.0000000434 * T**5
        )

    Nuta = nutation(T)
    ObliqAvg = obliqAvg(T) * S2R
    Obliq = ObliqAvg + Nuta["NutaObliq"]
    N = np.dot(np.dot(r1(-Obliq), r3(-Nuta["NutaEclp"])), r1(ObliqAvg))

    return N, Obliq


# 與js版驗算無誤
# print(nutaMx(1))
