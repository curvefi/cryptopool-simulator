#!/usr/bin/env python3

import json
import sys
import numpy as np
# from math import log
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = 'results.json'

with open(fname) as f:
    results = json.load(f)

As = set()
gammas = set()

for row in results['configuration']:
    As.add(row['A'])
    gammas.add(row['gamma'])

As = sorted(list(As))
gammas = sorted(list(gammas))

Z = np.zeros((len(gammas), len(As)))

for row in results['configuration']:
    # APY
    # liq_density
    # volume
    Z[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['slippage']
    # Z[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['APY']

fig, ax = plt.subplots()
plt.xscale('log')
plt.yscale('log')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
cbar = fig.colorbar(im, ax=ax)

ax.set_xlabel("A")
ax.set_ylabel("gamma")

ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200])
ax.set_xticklabels([10, 20, None, None, 50, None, None, None, None, 100, 200])

ax.set_yticks([2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
ax.set_yticklabels(["0.0002", "0.0005", "0.001", "0.002", "0.005", "0.01", "0.02"])

cbar.set_label("Slipage", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
