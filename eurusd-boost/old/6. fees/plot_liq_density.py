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
    As.add(row['mid_fee'])
    gammas.add(row['out_fee'])

As = sorted(list(As))
gammas = sorted(list(gammas))

Z = np.zeros((len(gammas), len(As)))

for row in results['configuration']:
    # APY
    # liq_density
    # volume
    Z[gammas.index(row['out_fee']), As.index(row['mid_fee'])] = row['Result']['liq_density']
    # Z[gammas.index(row['out_fee']), As.index(row['mid_fee'])] = row['Result']['APY']

fig, ax = plt.subplots()
plt.yscale('symlog')
plt.xscale('symlog')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
fig.colorbar(im, ax=ax)

plt.show()
