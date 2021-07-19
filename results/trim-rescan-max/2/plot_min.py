import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import log

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
    Z[gammas.index(row['gamma']), As.index(row['A'])] = log(row['Result']['slippage']-1.3)

print(Z.min(), Z.max())
fig, ax = plt.subplots()
plt.yscale('symlog')
plt.xscale('symlog')
im = ax.pcolormesh(As, gammas, Z, cmap=plt.get_cmap('jet'))
fig.colorbar(im, ax=ax)

plt.show()
