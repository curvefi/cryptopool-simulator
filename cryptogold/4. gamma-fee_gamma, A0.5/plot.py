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

x_axis = 'gamma'
y_axis = 'fee_gamma'

As = set()
gammas = set()

for row in results['configuration']:
    As.add(row[x_axis])
    gammas.add(row[y_axis])

As = sorted(list(As))
gammas = sorted(list(gammas))

Z = np.zeros((len(gammas), len(As)))

for row in results['configuration']:
    # APY
    # liq_density
    # volume
    Z[gammas.index(row[y_axis]), As.index(row[x_axis])] = row['Result']['slippage']

fig, ax = plt.subplots()
plt.xscale('log')
plt.yscale('log')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
cbar = fig.colorbar(im, ax=ax)

ax.set_xlabel(x_axis)
ax.set_ylabel(y_axis)

cbar.set_label("Slipage", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
