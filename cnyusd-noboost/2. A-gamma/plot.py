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
im.set_edgecolor('face')
cbar = fig.colorbar(im, ax=ax)

ax.set_xlabel("A")
ax.set_ylabel("gamma")

cbar.set_label("Price impact", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
