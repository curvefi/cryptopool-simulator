#!/usr/bin/env python3

import json
import sys
import pylab

X = "A"

if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = 'results.json'

with open(fname) as f:
    results = json.load(f)

x = []
Z = []

for row in results['configuration']:
    x.append(row[X])
    Z.append(row['Result']['APY'])

pylab.semilogx(x, Z)
pylab.xlabel("A")
pylab.ylabel("APY")
pylab.tight_layout()
pylab.show()
