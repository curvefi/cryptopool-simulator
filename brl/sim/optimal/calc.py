#!/usr/bin/env python3
from simulation import Curve

A = 230

# POOL sizes in USD
deposits = [
    5 * 10**6,
    10 * 10**6,
    20 * 10**6,
    50 * 10**6
]

# TX sizes in BRL
tx_sizes = [
    50_000,
    250_000,
    500_000,
    1 * 10**6,
    2_500_000,
    5 * 10**6
]

price = 0.18

for deposit in deposits:
    curve = Curve(A=A, D=(deposit * 10**18), n=2)
    print("Deposit:", deposit, "USD")
    print("---------------------")
    for size in tx_sizes:
        size_wad = int(price * size * 10**18)
        p0 = 10**18 / curve.dy(1, 0, 10**18)
        p1 = curve.dy(1, 0, size_wad) / size_wad
        impact = (p0 - p1) / p0 * 100
        print(f"Transaction size: {size} BRL  ->  Price impact: {impact:.4f}%")
    print()
