import json
import sys
import pylab

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = 'results.json'

with open(fname) as f:
    results = json.load(f)

x = []
Z = []

for row in results['configuration']:
    # APY
    # liq_density
    x.append(row['gamma'])
    Z.append(row['Result']['liq_density'])
    # Z.append(row['Result']['slippage'])

pylab.semilogx(x, Z)
pylab.show()
