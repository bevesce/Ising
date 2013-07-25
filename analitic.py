import matplotlib.pyplot as plt

from function import Function, FunctionBundle
import numpy as np
from frange import dl, dr


def M(H, B):
    J = 1
    K = J * B
    h = H * B
    m = (np.exp(K) * np.sinh(h)) / \
        ((np.exp(2 * K) * np.sinh(h) ** 2 + np.exp(-2 * K)) ** 0.5)
    return m

msl = []
for B in dl:
    mfs = []
    for H in dr:
        mfs.append(M(H, B))
    msl.append(Function(xs=dr, ys=mfs, xlabel='H', ylabel='M', label=str(B)))

ms = FunctionBundle(msl, dl)
ms.plot()
plt.figure()
