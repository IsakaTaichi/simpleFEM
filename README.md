# simpleFEM

四角要素の組み合わせで作る、構造物のFEM解析プログラムです。

## 必要なライブラリのインポート

import numpy as np
np.set_printoptions(threshold=10000)
# 要素描画用
import matplotlib.pyplot as plt
import matplotlib.patches as pat

# FEM計算
import simpleFEM
