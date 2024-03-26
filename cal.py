import numpy as np
from math import asin, atan2, cos, floor, sin, sqrt
from calpos import calPos
from numba import jit


pi = 3.141592653589793238462643
pi2 = 6.283185307179586476925287
S2R = 4.848136811095359935899141e-6


@jit
def calNewm(Y, Range):
    # 取癸卯曆1999年12月22日平冬至時間儒略日
    EpoSolsJd = 2451534.749
    # 採用癸卯曆首朔應，即十二月平朔距冬至的時間。與時憲曆用冬至次日夜半，我直接用冬至
    ChouConst = 15.68
    CloseOriginAd = 2000
    Solar = 365.2422
    Lunar = 29.530588853
    TermLeng = Solar / 12
    isNewm = True
    OriginAccum = (Y - CloseOriginAd) * Solar
    AvgSolsJd = EpoSolsJd + OriginAccum  # 歲前冬至
    AvgChouSd = (Lunar - OriginAccum % Lunar + ChouConst) % Lunar  # 首朔
    AcrJd = np.zeros(Range)
    for i in range(Range):
        AvgSd = AvgChouSd + (i + (0 if isNewm else 0.5)) * Lunar
        AvgJd = AvgSd + AvgSolsJd

        def delta(Jd):
            Sun = calPos("Sun", Jd)
            Moon = calPos("Moon", Jd)
            a = Sun[0] - Moon[0] - (0 if isNewm else pi)
            b = Moon[2] - Sun[2]
            if a < -7 / 4 * pi:
                a += pi2
            return a / b

        D = delta(AvgJd)
        AcrJd[i] = AvgJd
        while abs(D) > 1e-8:
            AcrJd[i] += D
            D = delta(AcrJd[i])
    np.savetxt("./newm.txt", AcrJd, fmt="%.8f", newline=", ")
    # return AcrJd


# print(calNewm(1024, 10))
