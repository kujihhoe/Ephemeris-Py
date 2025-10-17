# import inspect
# import sys
# import json
# from skyfield import almanac
# from skyfield.api import N, S, E, W, load, wgs84
from jplephem.spk import SPK
import numpy as np

from math import atan2, cos, sin, sqrt
from nutation import nutaMx
from precession import precessionMx

# print(inspect.getfile(SPK))
kernel = SPK.open("/Volumes/Pickman Ebooks/專書曆算/历表/DE/de440.bsp")
# # 使用[:, 1]来获取第二列的所有元素（索引从0开始，所以1表示第二列）
# column = matrix[:, 1]
# [2 5 8]
# 計算地金夾角：
# times=10 # 假設有10組位置
# a=np.zeros(len(times))
# for i in range (len(times)):
#     a[i]=np.arccos(ear[:,i].dot(ven[:,i])/np.linalg.norm(ear[:,i])/np.linalg.norm(ven[:,i]))*180/np.pi

# skyfield
# ts = load.timescale()
# eph = load('/Volumes/Pickman Ebooks/DE历表/de421.bsp')
# sun = eph['Sun']
# Beijing = wgs84.latlon(39.9062* N, 116.4284 * E)
# observer = eph['Earth'] + Beijing
# from skyfield import almanac_east_asia as almanac_ea

# t0 = ts.utc(2019, 12, 1)
# t1 = ts.utc(2019, 12, 31)
# t, tm = almanac.find_discrete(t0, t1, almanac_ea.solar_terms(eph))

# for tmi, ti in zip(tm, t):
#     print(tmi, almanac_ea.SOLAR_TERMS_ZHT[tmi], ti.utc_iso(' '))
#     # https://rhodesmill.org/skyfield/almanac.html#solar-terms
#     # The result t will be an array of times, and y will be integers in the range 0–23 which are each the index of a solar term. Localized Names for the solar terms in different East Asia languages are provided as SOLAR_TERMS_JP for Japanese, SOLAR_TERMS_VN for VietNamese, SOLAR_TERMS_ZHT for Traditional Chinese, and (as shown above) SOLAR_TERMS_ZHS for Simplified Chinese.


# print(kernel)
# DE440, 441:
# File type DAF/SPK and format LTL-IEEE with 14 segments:
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Mercury Barycenter (1)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Venus Barycenter (2)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Earth Barycenter (3)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Mars Barycenter (4)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Jupiter Barycenter (5)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Saturn Barycenter (6)
# 2287184.50..2688976.50  Type 2  Solar System Barycenter (0) -> Sun (10)
# 2287184.50..2688976.50  Type 2  Earth Barycenter (3) -> Moon (301)
# 2287184.50..2688976.50  Type 2  Earth Barycenter (3) -> Earth (399)
# 2287184.50..2688976.50  Type 2  Mercury Barycenter (1) -> Mercury (199)
# 2287184.50..2688976.50  Type 2  Venus Barycenter (2) -> Venus (299)
# 其中，
# Type 2 — positions stored as Chebyshev polynomials, with velocity derived by computing their derivative.
# Type 3 — positions and velocities both stored explicitly as Chebyshev polynomials.
# DE441的儒略日範圍：
# 441_1: -3100015.5 through 2440432.5
# 441_2: 2440400.5 through 8000016.5
# 440: 2287184.5 through 2688976.5

pi2 = 6.283185307179586476925287
S2R = 4.848136811095359935899141e-6
MMS2R = 4.848136811095359935899141e-13
c = 299792.458


def cal(a, b, jd):
    position, velocity = kernel[a, b].compute_and_differentiate(jd)
    return {"X": position, "V": velocity}


# print(cal(0, 10, 242111))

# def convert_to_array(numbers_str):
#     # 将字符串分割并转换为整数列表
#     return

# 用來從命令行獲取命令：
# if __Name__ == "__main__":
#     # 检查是否有足够的命令行参数
#     if len(sys.argv) > 1:
#         # 获取所有命令行参数（除了脚本名称），并将它们合并成一个字符串
#         input_str = ' '.join(sys.argv[1:])
#         # 将合并后的字符串转换为整数列表
#         numbers = [int(num) for num in input_str.split()]
#         # 打印结果
#         result = cal(numbers)
#         # 打印结果
#         print(result)
#     else:
#         print("Error: No input provided", file=sys.stderr)
# print(cal(0,10,2344211))
# 輸出的位置[x,y,z]爲ICRS，單位都是km。J2000地球平赤道面为平面，J2000平春分点方向为方向的直角坐标系。


B = np.array(
    [
        [0.99999999999999425, -7.078279744e-8, 8.05614894e-8],
        [7.078279478e-8, 0.99999999999999695, 3.306041454e-8],
        [-8.056149173e-8, -3.306040884e-8, 0.999999999999996208],
    ]
)


# 與VSOP驗算無誤。VSOP的x軸與DE完全一致，只差了幾十公里，但是y、z軸都不一樣，就是黃道赤道的區別。


def calXV(Name, jd):
    Planets = {
        "Sun": 10,
        "Mars": 4,
        "Jupiter": 5,
        "Saturn": 6,
        "Mercury": 199,
        "MercuryBary": 1,
        "Venus": 299,
        "VenusBary": 2,
    }
    if Name == "Sun" or Name == "Mars" or Name == "Jupiter" or Name == "Saturn":
        SB_S = cal(0, Planets[Name], jd)
        SB_EB = cal(0, 3, jd)
        EB_E = cal(3, 399, jd)
        X = np.subtract(SB_S["X"], np.add(SB_EB["X"], EB_E["X"]))
        V = np.subtract(SB_S["V"], np.add(SB_EB["V"], EB_E["V"]))
    elif Name == "Moon":
        EB_E = cal(3, 399, jd)
        EB_M = cal(3, 301, jd)
        X = np.subtract(EB_M["X"], EB_E["X"])
        V = np.subtract(EB_M["V"], EB_E["V"])
    elif Name == "Mercury" or Name == "Venus":
        SB_MB = cal(0, Planets[Name + "Bary"], jd)
        MB_M = cal(Planets[Name + "Bary"], Planets[Name], jd)
        SB_EB = cal(0, 3, jd)
        EB_E = cal(3, 399, jd)
        X = np.subtract(np.add(SB_MB["X"], MB_M["X"]), np.add(SB_EB["X"], EB_E["X"]))
        V = np.subtract(np.add(SB_MB["V"], MB_M["V"]), np.add(SB_EB["V"], EB_E["V"]))
    return X, V


def calPos(Name, Jd):
    def x2LonLat(X):
        Lon = atan2(X[1], X[0])
        Lat = atan2(X[2], sqrt(X[0] ** 2 + X[1] ** 2))
        return Lon, Lat

    def vec2dist(arr):
        return sqrt(arr[0] ** 2 + arr[1] ** 2 + arr[2] ** 2)

    def calLightAber(Jd, X):
        return Jd - vec2dist(X) / c / 86400

    def r1(a):
        return np.array([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]])

    T = (Jd - 2451545) / 36525  # 儒略世纪
    Xreal = calXV(Name, Jd)
    Jdr = calLightAber(Jd, Xreal[0])  # 推迟时
    X, V = calXV(Name, Jdr)  # 视位置
    X2000 = np.dot(B, X)
    V2000 = np.dot(B, V)
    # LonRaw calculation would be here, but is commented out in the provided code
    N, Obliq = nutaMx(T)
    P = precessionMx(T)  # 岁差矩阵 P(t)
    NP = np.dot(N, P)
    R1Eps = r1(Obliq)
    Equa = np.dot(NP, X2000)
    Equa1 = np.dot(NP, V2000)
    Eclp = np.dot(R1Eps, Equa)
    Eclp1 = np.dot(R1Eps, Equa1)
    Lon1 = (Eclp[0] * Eclp1[1] - Eclp[1] * Eclp1[0]) / (Eclp[0] ** 2 + Eclp[1] ** 2)
    # 转换
    EquaLon, EquaLat = x2LonLat(Equa)
    EquaLon = (EquaLon + pi2) % pi2
    Lon, Lat = x2LonLat(Eclp)
    Lon = (Lon + pi2) % pi2

    return Lon, Lat, Lon1, EquaLon, EquaLat


# print(calPos("Sun", 1433133))
