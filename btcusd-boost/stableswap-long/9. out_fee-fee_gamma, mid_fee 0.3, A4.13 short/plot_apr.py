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

x_axis = 'out_fee'
y_axis = 'fee_gamma'

for row in results['configuration']:
    As.add(row[x_axis])
    gammas.add(row[y_axis])

As = sorted(list(As))
gammas = sorted(list(gammas))

Z = np.zeros((len(gammas), len(As)))

for row in results['configuration']:
    apy = row['Result']['APY']
    apr = ((1 + apy)**(1 / (365 * 86400)) - 1) * 365 * 86400
    Z[gammas.index(row[y_axis]), As.index(row[x_axis])] = apr - row['boost_rate']

fig, ax = plt.subplots()
plt.yscale('log')
plt.xscale('log')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
im.set_edgecolor('face')
cbar = fig.colorbar(im, ax=ax)

ax.set_xlabel("out_fee")
ax.set_ylabel("fee_gamma")
cbar.set_label("APY", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
