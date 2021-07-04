import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import log

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = 'results.json'

with open(fname) as f:
    results = json.load(f)

mid_fee = set()
out_fee = set()

for row in results['configuration']:
    mid_fee.add(row['fee_gamma'])
    out_fee.add(row['out_fee'])

mid_fee = sorted(list(mid_fee))
out_fee = sorted(list(out_fee))

Z = np.zeros((len(mid_fee), len(out_fee)))

for row in results['configuration']:
    # APY
    # liq_density
    # volume
    print(row['Result']['slippage'])
    Z[mid_fee.index(row['fee_gamma']), out_fee.index(row['out_fee'])] = log(row['Result']['slippage'])

fig, ax = plt.subplots()
plt.yscale('symlog')
plt.xscale('symlog')
im = ax.pcolormesh(out_fee, mid_fee, Z, cmap=plt.get_cmap('jet'))
fig.colorbar(im, ax=ax)

plt.show()
