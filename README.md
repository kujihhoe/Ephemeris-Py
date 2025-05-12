# Ephemeris-Py
DE历表计算朔闰节气

|          | 模型                                                         | 代碼來源                                                     |
| -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 總體方法 | 廖育棟                                                       | 《[月相和二十四節氣的計算](https://ytliu0.github.io/ChineseCalendar/computation_simp.html)》、*[Calculations in Star Charts](https://ytliu0.github.io/starCharts/)* |
| 曆表     | DE441、440                                                   | [jplephem](https://pypi.org/project/jplephem/)               |
| 歲差     | Vondrak 等（2011）                                           | [Vondrak](https://github.com/dreamalligator/vondrak)         |
| 章動     | [IAU2000A](http://asa.usno.navy.mil/SecM/Glossary.html#nutation) | [python-novas](https://github.com/brandon-rhodes/python-novas)/Cdist/nutation.c，由 GPT 轉換爲 python |

先安裝jplephem：

```
pip install jplephem
```

