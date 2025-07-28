#!/usr/bin/env python3
import pylab
from simulation import Curve
import numpy as np

A = 100
deposit = 1000 * 10**18

curve = Curve(A=A, D=deposit, n=2)
xx = curve.x[:]

prices = []
impacts = []

for swap_size in np.linspace(1e-5, 0.5, 300):
    curve.x = xx[:]
    dx = int(swap_size * curve.x[0])
    curve.exchange(0, 1, dx)
    dx = int(min(curve.x) / 1e8)
    dy = curve.dy(0, 1, dx)
    p_1 = dx / dy

    dx_impact = dx * 100

    curve.exchange(0, 1, dx_impact)
    dy = curve.dy(0, 1, dx)
    p_2 = dx / dy

    impact = (p_2 - p_1) / (dx / 1e18)

    prices.append(p_1)
    impacts.append(impact)

dp = np.array(prices) - 1
pylab.plot(dp, impacts, c='black')
pylab.plot(dp, 1 / (2 * A) * (1 + (dp * A)**2), c='red')
pylab.show()
