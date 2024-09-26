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

x_axis = 'mid_fee'
y_axis = 'out_fee'

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
    Z[gammas.index(row[y_axis]), As.index(row[x_axis])] = row['Result']['liq_density']

fig, ax = plt.subplots()
plt.yscale('log')
plt.xscale('log')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
fig.colorbar(im, ax=ax)

plt.show()
