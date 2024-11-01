#!/usr/bin/env python3

import json
import sys
import numpy as np
# from math import log
import pylab

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

slippage = np.zeros((len(gammas), len(As)))
APY = np.zeros((len(gammas), len(As)))
liq_density = np.zeros((len(gammas), len(As)))

for row in results['configuration']:
    # APY
    # liq_density
    # volume
    slippage[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['slippage']
    APY[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['APY']
    liq_density[gammas.index(row['gamma']), As.index(row['A'])] = row['Result']['liq_density']

gamma_ix = slippage.argmin(axis=0)
result = [APY[g, a] for (g, a) in zip(gamma_ix, range(len(As)))]

pylab.semilogx(As, result)
pylab.show()
