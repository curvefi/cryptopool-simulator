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
    Z[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['volume']
    # Z[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['APY']

fig, ax = plt.subplots()
plt.yscale('log')
plt.xscale('log')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
im.set_edgecolor('face')
cbar = fig.colorbar(im, ax=ax)

ax.set_xlabel("A")
ax.set_ylabel("gamma")
cbar.set_label("Trading volume", rotation=270, labelpad=15)

ax.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40])
ax.set_xticklabels([0.2, None, None, 0.5, None, None, None, None, 1, 2, None, None, 5, None, None, None, None, 10, 20,
                    None, None])

plt.tight_layout()
plt.show()
