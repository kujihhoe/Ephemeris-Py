import numpy as np
from calpos import calPos


pi = 3.141592653589793238462643
pi2 = 6.283185307179586476925287
S2R = 4.848136811095359935899141e-6
D2R = 0.0174532925199432957692369
# 取癸卯曆1999年12月22日平冬至時間儒略日
EpoSolsJd = 2451534.749
# 採用癸卯曆首朔應，即十二月平朔距冬至的時間。與時憲曆用冬至次日夜半，我直接用冬至
ChouConst = 15.68
CloseOriginAd = 2000
Solar = 365.2422
Lunar = 29.530588853
TermLeng = Solar / 12


# isNewm True：算朔，False：望
# isTerm True：算中氣，False：算節氣
def calNewm(Y, Range, isNewm):
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
            if a > .8 * pi:
                a -= pi2
            elif a < -.8 * pi:
                a += pi2
            return a / b
# syzygy：只有a<-pi的情况
# S 6    6     3    3.2   0.2  0.2   
# M 3    3.2    6    6     3    3.2
# a 3    2.8   -3   -2.8  -2.8  -3
#   -pi
# r -.1 -.3    -6.1 -5.9  
        D = delta(AvgJd)
        AcrJd[i] = AvgJd
        while abs(D) > 1e-8:
            AcrJd[i] += D
            D = delta(AcrJd[i])
    np.savetxt("../newm_2.txt", AcrJd, fmt="%.8f", newline=", ")
    # return AcrJd

def calTerm(Y, Range, isTerm):
    OriginAccum = (Y - CloseOriginAd) * Solar
    AvgSolsJd = EpoSolsJd + OriginAccum  # 歲前冬至
    AcrTermJd = np.zeros(Range)
    for i in range(Range):
        TermLon = D2R * (((2 * i + (2 if isTerm else 1)) * 15 + 270) % 360)
        AvgTermSd = (i + 1) * TermLeng
        AvgTermJd = AvgTermSd + AvgSolsJd

        def delta(Jd):
            Sun = calPos("Sun", Jd)
            a = TermLon - Sun[0]
            if a < -7 / 4 * pi:
                a += pi2
            b = Sun[2]
            return a / b

        D = delta(AvgTermJd)
        AcrTermJd[i] = AvgTermJd
        while abs(D) > 1e-8:
            AcrTermJd[i] += D
            D = delta(AcrTermJd[i])
    np.savetxt("term1_plus.txt", AcrTermJd, fmt="%.8f", newline=", ")
    # return AcrJd


print(calNewm(1600, 11140, True))
# print(calTerm(-676, 5, True))

