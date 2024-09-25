#!/usr/bin/env python3

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
Z2 = []

for row in results['configuration']:
    gammas.append(row['ma_half_time'])
    Z2.append(row['Result']['APY'])

# pylab.semilogx(gammas, Z1)
pylab.semilogx(gammas, Z2)
pylab.show()
