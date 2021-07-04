import json
import sys
import pylab

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = 'results.json'

with open(fname) as f:
    results = json.load(f)

gammas = []
Z1 = []
Z2 = []
Z3 = []
Z4 = []

for row in results['configuration']:
    # APY
    # liq_density
    gammas.append(row['gamma'])
    # Z1.append(row['Result']['liq_density'])
    Z2.append(row['Result']['slippage'])

# pylab.semilogx(gammas, Z1)
pylab.semilogx(gammas, Z2)
pylab.show()
